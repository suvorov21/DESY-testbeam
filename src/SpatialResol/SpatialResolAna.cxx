//#include <iostream>

#include "SpatialResolAna.hxx"

SpatialResolAna::SpatialResolAna(int argc, char** argv): AnalysisBase(argc, argv) {

  if (_iteration == -1) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna(). Iteration should be defined as a input param" << std::endl;
    exit(1);
  }

}

bool SpatialResolAna::Initialize() {
  std::cout << "*****************************************" << std::endl;
  std::cout << "***   Spatial resolution analysis    ****" << std::endl;
  std::cout << "*****************************************" << std::endl;

  std::cout << "Batch mode    :   " << _batch         << std::endl;
  std::cout << "Verbosity     :   " << _verbose       << std::endl;
  std::cout << "Iteration     :   " << _iteration     << std::endl;


  AnalysisBase::Initialize();

  std::cout << "Initializing spatial resolution ana......";

  gErrorIgnoreLevel = kSysError;

  if (_iteration) {
    // TODO Read the PRF params
    TString prev_file_name = _file_out_name;
    if (_iteration < 11) {
      prev_file_name = prev_file_name(0, prev_file_name.Last('.') - 1);
      prev_file_name +=  + TString::Itoa(_iteration - 1, 10);
      prev_file_name += ".root";

      _Prev_iter_file = new TFile(prev_file_name.Data(), "READ");
      _PRF_function   = (TF1*)_Prev_iter_file->Get("PRF_function");
      auto uncertainty_hist = (TH1F*)_Prev_iter_file->Get("resol_final");
      _uncertainty = 0;
      for (auto i = 2; i < geom::nPadx; ++i)
        _uncertainty += uncertainty_hist->GetBinContent(i) / (geom::nPadx - 2);

      Int_t read_var;
      auto event_tree = (TTree*)_Prev_iter_file->Get("EventTree");
      event_tree->SetBranchAddress("PassedEvents",    &read_var);
      std::vector<Int_t> vec;
      vec.clear();
      for (auto i = 0; i < event_tree->GetEntries(); ++i) {
        event_tree->GetEntry(i);
        vec.push_back(read_var);
      }
      this->SetEventList(vec);

    }
    if (!_PRF_function) {
      std::cerr << "ERROR. SpatialResolAna::Initialize(). PRF function is not specified" << std::endl;
      std::cerr << "Search in " << prev_file_name << std::endl;
      exit(1);
    }
  } else {
    /*_PRF_function = new TF1("PRF_function",
      " [0] * exp(-4*(1-[1])*TMath::Power(x/[2], 2.)) / (1+4 * [1] * TMath::Power(x/[2], 2.) )", prf_min, prf_max);
    _PRF_function->SetParName(0, "Const");
    _PRF_function->SetParName(1, "r");
    _PRF_function->SetParName(2, "w");

    auto c = 0.8;
    auto r = 0.5;
    auto s = 0.7;
    _PRF_function->SetParameters(c, r, s);*/

    _PRF_function  = new TF1("PRF_function", "[0]*(1+[1]*x*x + [2] * x*x*x*x) / (1+[3]*x*x+[4]*x*x*x*x)", prf_min, prf_max);
    _PRF_function->SetParName(0, "Const");
    _PRF_function->SetParName(1, "a2");
    _PRF_function->SetParName(2, "a4");
    _PRF_function->SetParName(3, "b2");
    _PRF_function->SetParName(4, "b4");

    double co = 0.83;
    double a2 = -0.84;
    double a4 = 0.83;
    double b2 = 25.2;
    double b4 = 0.83;
    _PRF_function->SetParameters(co, a2, a4, b2, b4);
  }

  // Initialise histoes and graphs
  _PRF_histo = new TH2F("PRF_histo","", prf_bin, prf_min, prf_max, 102,0.,1.02);
  _PRF_histo_2pad = new TH2F("PRF_histo_2pad","", prf_bin, prf_min, prf_max, 102,0.,1.02);
  _PRF_histo_3pad = new TH2F("PRF_histo_3pad","", prf_bin, prf_min, prf_max, 102,0.,1.02);
  _PRF_histo_4pad = new TH2F("PRF_histo_4pad","", prf_bin, prf_min, prf_max, 102,0.,1.02);

  _PRF_graph = new TGraphErrors();
  _PRF_graph->SetName("PRF_graph");

  for (auto j = 0; j < geom::nPadx; ++j) {
    _resol_col_hist[j]  = new TH1F(Form("resol_histo_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_except[j]  = new TH1F(Form("resol_histo1__%i", j), "", resol_bin, resol_min, resol_max);

    _resol_col_hist_2pad[j]  = new TH1F(Form("resol_histo_2pad_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_2pad_except[j]  = new TH1F(Form("resol_histo1_2pad_%i", j), "", resol_bin, resol_min, resol_max);

    _resol_col_hist_3pad[j]  = new TH1F(Form("resol_histo_3pad_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_3pad_except[j]  = new TH1F(Form("resol_histo1_3pad_%i", j), "", resol_bin, resol_min, resol_max);
  }

  _residual_sigma_unbiased  = new TH1F("resol", "", geom::nPadx, 0., geom::nPadx);
  _residual_sigma_biased    = new TH1F("resol1", "", geom::nPadx, 0., geom::nPadx);
  _residual_sigma           = new TH1F("resol_final", "", geom::nPadx, 0., geom::nPadx);

  _residual_mean            = new TH1F("mean", "", geom::nPadx, 0., geom::nPadx);

  _Chi2_track = new TH1F("Chi2_Track", "", 100, 0., 3.);

  // schedule the output for writing
  _output_vector.push_back(_PRF_function);
  _output_vector.push_back(_PRF_histo);
  _output_vector.push_back(_PRF_histo_2pad);
  _output_vector.push_back(_PRF_histo_3pad);
  _output_vector.push_back(_PRF_histo_4pad);
  _output_vector.push_back(_PRF_graph);

  _output_vector.push_back(_residual_sigma);
  _output_vector.push_back(_residual_mean);

  _output_vector.push_back(_residual_sigma_biased);
  _output_vector.push_back(_residual_sigma_unbiased);

  _output_vector.push_back(_Chi2_track);

  for (auto j = 0; j < geom::nPadx; ++j) {
    _output_vector.push_back(_resol_col_hist[j]);
  }

  _passed_events.clear();

  std::cout << "done" << std::endl;
  std::cout << "      PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula() << "  with ";
  for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
    std::cout << "  " << _PRF_function->GetParameter(i) << ",";
  std::cout << std::endl;

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  return true;
}

bool SpatialResolAna::ProcessEvent(const TEvent* event) {

  bool passed = false;

  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    if(sel::GetNonZeroCols(event  , trackId).size() != geom::nPadx) return false;
    if(sel::GetNonZeroRows(event  , trackId).size()>5) return false;
    if(sel::GetFitQuality(event   , trackId)>1.0e6) return false;

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    passed = true;

    if (_verbose == 2)
      std::cout << "Track id = " << trackId << std::endl;

    TGraphErrors* track_gr = new TGraphErrors();
    TGraphErrors* track_m = new TGraphErrors();

    TGraphErrors* track_1[36];
    for (int i = 0; i < 36; ++i) {
      track_1[i] = new TGraphErrors();
    }

    int cluster[36];
    int cluster_N[36];
    double true_track[36];
    double cluster_mean[36];
    float charge_max[36];

    // At the moment ommit first and last column
    // first loop over column
    for (uint colID = 0; colID < track->GetCols().size(); ++colID) {
      // FIXIT check this call
      auto it_x = track->GetColHits(colID)[0]->GetCol();
      cluster[it_x]     = 0;
      cluster_N[it_x]   = 0;
      charge_max[it_x]  = 0;
      true_track[it_x]  = -999.;
      cluster_mean[it_x] = 0.;

      TH1F* cluster_h = new TH1F("cluster", "", geom::nPady, -1.*geom::MM_dy - geom::dy, geom::MM_dy + geom::dy);

      for (uint rowID = 0; rowID < track->GetColHits(it_x).size(); ++rowID) {
        if (!track->GetColHits(it_x)[rowID])
          continue;

        auto it_y = track->GetColHits(it_x)[rowID]->GetRow();
        auto q = track->GetColHits(colID)[rowID]->GetQ();

        cluster[it_x] += q;
        ++cluster_N[it_x];
        cluster_h->Fill(geom::y_pos[it_y], q);
        if (charge_max[it_x] < q)
          charge_max[it_x] = q;
      } // end of loop over rows

      cluster_mean[it_x] = cluster_h->GetMean();

      delete cluster_h;

      if (!cluster[it_x] || !charge_max[it_x])
        continue;

      // minimise chi2 to extimate true_track
      double chi2_min     = 1e9;
      double scan_y       = cluster_mean[it_x] - scan_delta;

      if (_iteration > 0 && cluster_N[it_x] != 1) {
        for (Int_t scanId = 0; scanId < scan_Nsteps; ++scanId) {
          double chi2 = 0;
          for (uint rowID = 0; rowID < track->GetColHits(it_x).size(); ++rowID) {
            auto q      = track->GetColHits(colID)[rowID]->GetQ();
            auto it_y   = track->GetColHits(colID)[rowID]->GetRow();
            if (!q)
              continue;

            double a = 1. * q / cluster[it_x];
            double center_pad_y = geom::y_pos[it_y];
            double part = (a - _PRF_function->Eval(scan_y - center_pad_y));
            double c = 1.*cluster[it_x];
            double b = 1.*q;
            part *= c*c;
            part /= c*sqrt(b) + b*sqrt(c);
            part *= part;

            chi2 += part;
          }

          if (chi2 < chi2_min) {
            chi2_min = chi2;
            true_track[it_x] = scan_y;
          }

          scan_y += scan_step;
        }
      } else
      true_track[it_x] = cluster_mean[it_x];

      double x = geom::x_pos[it_x];

      track_gr->SetPoint(track_gr->GetN(), x, true_track[it_x]);
      for (int i = 1; i < geom::nPadx - 1; ++i) {
        if (i != it_x)
          track_1[i]->SetPoint(track_1[i]->GetN(), x, true_track[it_x]);
      }

      track_m->SetPoint(track_m->GetN(), x, cluster_mean[it_x]);
      double error;
      if (cluster_N[it_x] == 1)
        error = one_pad_error;
      else {
        if (_iteration == 0)
          error = default_error;
        else
          error = _uncertainty;
      }
      track_gr->SetPointError(track_gr->GetN()-1, 0., error);
      for (int i = 1; i < geom::nPadx; ++i) {
        if (i != it_x)
          track_1[i]->SetPointError(track_1[i]->GetN() - 1, 0., error);
      }
    } // loop over i

    track_gr->Fit("pol1", "Q");
    TF1* fit = track_gr->GetFunction("pol1");

        //track_m->Fit("pol1", "Q");
        //TF1* fit1 = track_m->GetFunction("pol1");

    if (!fit)
      continue;

    double quality = fit->GetChisquare() / fit->GetNDF();
    double k = fit->GetParameter(1);
    double b = fit->GetParameter(0);

    _Chi2_track->Fill(quality);

    double k1[36];
    double b1[36];

    for (int i = 1; i < geom::nPadx; ++i) {
      track_1[i]->Fit("pol1", "Q");
      TF1* fit1 = track_1[i]->GetFunction("pol1");
      if (!fit1)
        continue;
      k1[i] = fit1->GetParameter(1);
      b1[i] = fit1->GetParameter(0);
    }

    // second loop over columns
    for (uint colID = 0; colID < track->GetCols().size(); ++colID) {
      // FIXIT check this call
      auto it_x = track->GetColHits(colID)[0]->GetCol();

      if (true_track[it_x]  == -999.)
        continue;

      double x    = geom::x_pos[it_x];
      double track_fit_y  = k * x + b;
      double track_fit_y1 = k1[it_x] * x + b1[it_x];

      // fill SR
      _resol_col_hist[it_x]->Fill(true_track[it_x] - track_fit_y);
      _resol_col_hist_except[it_x]->Fill(true_track[it_x] - track_fit_y1);

      if (cluster_N[it_x] == 2) {
        _resol_col_hist[it_x]->Fill(true_track[it_x] - track_fit_y);
        _resol_col_hist_2pad_except[it_x]->Fill(true_track[it_x] - track_fit_y1);
      } else if (cluster_N[it_x] == 3) {
        _resol_col_hist_3pad[it_x]->Fill(true_track[it_x] - track_fit_y);
        _resol_col_hist_3pad_except[it_x]->Fill(true_track[it_x] - track_fit_y1);
      }

      if (cluster_N[it_x] == 1)
        continue;
      // Fill PRF
      for (uint rowID = 0; rowID < track->GetColHits(it_x).size(); ++rowID) {
        if (!track->GetColHits(it_x)[rowID])
          continue;

        auto it_y = track->GetColHits(it_x)[rowID]->GetRow();
        auto q = track->GetColHits(colID)[rowID]->GetQ();

        if (!cluster[it_x] || !q)
          continue;

        double charge = 1. * q / cluster[it_x];
        double center_pad_y = geom::y_pos[it_y];

            // fill PRF
        _PRF_histo->Fill(center_pad_y - track_fit_y, charge);

        if (cluster_N[it_x] == 2)
          _PRF_histo_2pad->Fill(center_pad_y - track_fit_y, charge);
        else if (cluster_N[it_x] == 3)
          _PRF_histo_3pad->Fill(center_pad_y - track_fit_y, charge);
        else if (cluster_N[it_x] == 4)
          _PRF_histo_4pad->Fill(center_pad_y - track_fit_y, charge);
      }
    } // loop over colums
    delete track_gr;
    delete track_m;
    for (int i = 0; i < 36; ++i) {
      delete track_1[i];
    }
  } // loop over tracks

  if (passed)
    _passed_events.push_back(event->GetID());

  return true;
}

bool SpatialResolAna::WriteOutput() {
  std::cout << "Postprocessing histoes for writing.......";

  // Output histoes postprocession
  for (auto i = 1; i < _PRF_histo->GetXaxis()->GetNbins(); ++i) {

    TH1D* temp_h = _PRF_histo-> ProjectionY(Form("projections_bin_%i", i), i, i);

    double x = _PRF_histo->GetXaxis()->GetBinCenter(i);
    double y = temp_h->GetBinCenter(temp_h->GetMaximumBin());

    float start = -1.;
    float end   = -1.;
    float max = temp_h->GetMaximum();

    for (Int_t bin = 0; bin < temp_h->GetXaxis()->GetNbins(); ++bin) {
      if (start == -1. && temp_h->GetBinContent(bin) >= max / 2.)
        start = temp_h->GetBinCenter(bin);

      if (end == -1. && start != -1. && temp_h->GetBinContent(bin) <= max / 2.)
        end = temp_h->GetBinCenter(bin);
    }

    float e = end - start;

    _PRF_graph->SetPoint(_PRF_graph->GetN(), x, y);
    _PRF_graph->SetPointError(_PRF_graph->GetN()-1, 0, e/2.);
  } // end of PRF histo profiling

  for (auto i = 1; i < geom::nPadx - 1; ++i) {
    _resol_col_hist[i]->Fit("gaus", "Q");
    _resol_col_hist_2pad[i]->Fit("gaus", "Q");
    _resol_col_hist_3pad[i]->Fit("gaus", "Q");

    _resol_col_hist_except[i]->Fit("gaus", "Q");
    _resol_col_hist_2pad_except[i]->Fit("gaus", "Q");
    _resol_col_hist_3pad_except[i]->Fit("gaus", "Q");

    auto mean     = _resol_col_hist[i]->GetFunction("gaus")->GetParameter(1);
    auto sigma    = _resol_col_hist[i]->GetFunction("gaus")->GetParameter(2);
    auto sigma_ex = _resol_col_hist_except[i]->GetFunction("gaus")->GetParameter(2);

    _residual_sigma_biased->SetBinContent(i+1, sigma);
    _residual_sigma_unbiased->SetBinContent(i+1, sigma_ex);
    _residual_sigma->SetBinContent(i+1, sqrt(sigma * sigma_ex));

    _residual_mean->SetBinContent(i+1, mean);
  }

  _PRF_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);

  // Output histoes postprocession done

  std::cout << "done" << std::endl;
  std::cout << "      PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula() << "  with ";
  for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
    std::cout << "  " << _PRF_function->GetParameter(i) << ",";
  std::cout << "  Chi2/NDF " << _PRF_graph->GetFunction("PRF_function")->GetChisquare()
            << "/" << _PRF_graph->GetFunction("PRF_function")->GetNDF() << std::endl;
  std::cout << std::endl;

  // Write objects
  AnalysisBase::WriteOutput();

  if (!_file_out)
    return true;

  std::cout << "Writing spatial analisis output..........";

  auto file = new TFile(_file_out_name.Data(), "UPDATE");
  // write
  auto tree = new TTree("EventTree", "");
  Int_t var = 0;
  tree->Branch("PassedEvents",     &var);
  for (uint i = 0; i < _passed_events.size(); ++i) {
    var = _passed_events[i];
    tree->Fill();
  }
  tree->Write("", TObject::kOverwrite);
  file->Close();

  if (_Prev_iter_file)
    _Prev_iter_file->Close();

  std::cout << "done" << std::endl;
  return true;
}

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna(argc, argv);
  ana->Initialize();
  ana->Loop(ana->GetEventList());
  ana->WriteOutput();

  return 1;
}
