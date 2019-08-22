//#include <iostream>

#include "SpatialResolAna.hxx"

#include "Math/Functor.h"
#include "Fit/Fitter.h"

//! verbosity level
//! 1 (default) print progress, memory usage, time consumption
//! 2           print event and track number
//! 3           print fit details
//! 4           print PRF details

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

    /*double g  = 3. * 0.01;
    double delta  = g / 2.;
    double a      = 1.;
    double b      = 0.;


    double co = 1.;
    double a2 = -1.* TMath::Power(2/delta, 2) * (1+a);
    double a4 = TMath::Power(2./delta, 4) * a;
    double b2 = TMath::Power(2./g, 2) * (1-b-2*(1+a)*TMath::Power(g/delta, 2) + 2.*a*TMath::Power(g/delta, 4));
    double b4 = TMath::Power(2./g, 4)*b;*/
    double co = 1.;
    double a2 = 2.35167e3;
    double a4 = 6.78962e7;
    double b2 = 3.36748e3;
    double b4 = 6.45311e8;


    _PRF_function->SetParameters(co, a2, a4, b2, b4);
    _PRF_function->FixParameter(0, 1.);
  }

  _do_arc_fit        = true;
  _do_full_track_fit = false;

  _circle_function_up = new TF1("circle_up", "-sqrt([0]*[0] - TMath::Power(x+0.198 - [1] * [0], 2)) + [0] * sqrt(1-[1]*[1]) + [2]", -0.5, 0.5);
  _circle_function_up->SetParName(0, "radius");
  _circle_function_up->SetParName(1, "sin(alpha)");
  _circle_function_up->SetParName(2, "target");

  _circle_function_dn = new TF1("circle_dn", "sqrt([0]*[0] - TMath::Power(x+0.198 - [1] * [0], 2)) - [0] * sqrt(1-[1]*[1]) + [2]", -0.5, 0.5);
  _circle_function_dn->SetParName(0, "radius");
  _circle_function_dn->SetParName(1, "sin(alpha)");
  _circle_function_dn->SetParName(2, "target");

  _qulity_ratio   = new TH1F("quality_ratio", "Quality ratio arc/linear", 100, 0., 2.);
  _chi2_ratio     = new TH1F("chi_ratio", "Chi2 ratio arc/linear", 100, 0., 2.);

  _mom_reco   = new TH1F("mom_reco", "", 2000, -6., 6.);
  _pos_reco   = new TH1F("pos_reco", "", 8000, -0.2, 0.2);
  _ang_reco   = new TH1F("ang_reco", "", 3000, -0.3, 0.3);

  // Initialise histoes and graphs
  _PRF_histo = new TH2F("PRF_histo","", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _PRF_histo_2pad = new TH2F("PRF_histo_2pad","", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _PRF_histo_3pad = new TH2F("PRF_histo_3pad","", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _PRF_histo_4pad = new TH2F("PRF_histo_4pad","", prf_bin, prf_min, prf_max, 150,0.,1.5);

  for (auto i = 0; i < geom::nPadx; ++i)
    _PRF_histo_col[i] = new TH2F(Form("PRF_histo_col_%i", i),"", prf_bin, prf_min, prf_max, 150,0.,1.5);

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
  _output_vector.push_back(_qulity_ratio);
  _output_vector.push_back(_chi2_ratio);
  _output_vector.push_back(_mom_reco);
  _output_vector.push_back(_pos_reco);
  _output_vector.push_back(_ang_reco);

  for (auto j = 0; j < geom::nPadx; ++j) {
    _output_vector.push_back(_resol_col_hist[j]);
  }

  for (auto j = 0; j < geom::nPadx; ++j) {
    _output_vector.push_back(_PRF_histo_col[j]);
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

  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    if(sel::GetNonZeroCols(track).size() != geom::nPadx) return false;
    if(sel::GetNonZeroRows(track).size() > 8) return false;
    if(sel::GetColsMaxSep(track) > 6) return false;
    if(sel::GetColsMaxGap(track) > 1) return false;

    if(_test_mode) DrawSelection(event,trackId);

    _store_event = true;

    if (_verbose > 1)
      std::cout << "Track id = " << trackId << std::endl;

    TGraphErrors* track_gr = new TGraphErrors();
    TGraphErrors* track_m = new TGraphErrors();

    TGraphErrors* track_1[geom::nPadx];
    for (int i = 0; i < geom::nPadx; ++i) {
      track_1[i] = new TGraphErrors();
    }

    int cluster[geom::nPadx];
    int cluster_N[geom::nPadx];
    double track_pos[geom::nPadx];
    double cluster_mean[geom::nPadx];
    float charge_max[geom::nPadx];
    double a_peak[geom::nPadx];

    // At the moment ommit first and last column
    // first loop over columns
    for (auto row:track->GetCols()) {
      if (!row[0])
        continue;
      auto it_x = row[0]->GetCol();
      // exlude 1st/last column
      if (it_x == 0 || it_x == geom::nPadx-1)
        continue;

      cluster[it_x]       = 0;
      cluster_N[it_x]     = 0;
      charge_max[it_x]    = 0;
      track_pos[it_x]     = -999.;
      cluster_mean[it_x]  = 0.;
      a_peak[it_x]        = 0.;

      TH1F* cluster_h = new TH1F("cluster", "", geom::nPady, -1.*geom::MM_dy - geom::dy, geom::MM_dy + geom::dy);

      // loop over rows
      for (auto pad:row) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow();
        auto q = pad->GetQ();

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

      if (/*_iteration > 0 && */cluster_N[it_x] != 1) {
        for (Int_t scanId = 0; scanId < scan_Nsteps; ++scanId) {
          // TODO implement ROOT fitter instead of manual fit
          double chi2 = 0;

          double a_nom = 0.;
          double a_den = 0.;
          double a_tot = 0.;

          for (auto pad:row) {
            auto q      = pad->GetQ();
            auto it_y   = pad->GetRow();
            if (!q)
              continue;
            double center_pad_y = geom::y_pos[it_y];

            a_nom += _PRF_function->Eval(scan_y - center_pad_y);
            a_den += TMath::Power(_PRF_function->Eval(scan_y - center_pad_y), 2) / q;
          }
          a_tot = a_nom / a_den;

          for (auto pad:row) {
            auto q      = pad->GetQ();
            auto it_y   = pad->GetRow();
            if (!q)
              continue;
            double center_pad_y = geom::y_pos[it_y];

            double part = (q - a_tot*_PRF_function->Eval(scan_y - center_pad_y));
            part *= part;
            part /= q;

            chi2 += part;
          }

          if (chi2 < chi2_min) {
            chi2_min = chi2;
            track_pos[it_x] = scan_y;
            a_peak[it_x] = a_tot;
          }

          scan_y += scan_step;
        }
      } else
        track_pos[it_x] = cluster_mean[it_x];

      double x = geom::x_pos[it_x];

      track_gr->SetPoint(track_gr->GetN(), x, track_pos[it_x]);
      for (int i = 1; i < geom::nPadx - 1; ++i) {
        if (i != it_x)
          track_1[i]->SetPoint(track_1[i]->GetN(), x, track_pos[it_x]);
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
      for (int i = 1; i < geom::nPadx-1; ++i) {
        if (i != it_x)
          track_1[i]->SetPointError(track_1[i]->GetN() - 1, 0., error);
      }
    } // loop over columns

    TF1* fit;
    TF1* lin_fit;
    TString func = "";
    if (!_do_arc_fit) {
      func = "pol1";
      track_gr->Fit(func, "Q");
      fit = track_gr->GetFunction(func);
    } else {
/*
      auto chi2Function = [&](const Double_t *par) {
        //minimisation function computing the sum of squares of residuals
        // looping at the graph points
        Int_t np = track_gr->GetN();
        Double_t f = 0;
        Double_t *x = track_gr->GetX();
        Double_t *y = track_gr->GetY();
        Double_t *ye = track_gr->GetEY();
        for (Int_t i=0;i<np;i++) {
           Double_t u = x[i] - par[0] * par[1] + 0.198;
           Double_t v = y[i] - par[0]*sqrt(1+par[1]*par[1])+par[2];
           Double_t dr = par[0] - std::sqrt(u*u+v*v);
           Double_t re = ye[i] * (y[i] - par[0]/sqrt(1+par[1]*par[1]) + par[2]) / par[0];
           dr /= re;
           f += dr*dr;
        }
        return f;
      };

      ROOT::Math::Functor fcn(chi2Function,3);
      ROOT::Fit::Fitter  fitter;

      double pStart[3] = {80., 0., 0.};
      fitter.SetFCN(fcn, pStart);
      fitter.Config().ParSettings(0).SetName("R");
      fitter.Config().ParSettings(1).SetName("tan");
      fitter.Config().ParSettings(2).SetName("c");

      bool ok = fitter.FitFCN();
      const ROOT::Fit::FitResult & result = fitter.Result();
      result.Print(std::cout);
*/

      // *********** end of debug
      _circle_function_dn->SetParameters(80., 0., 0.);
      //_circle_function_dn->SetParLimits(0, 5., 1000.);
      _circle_function_dn->SetParLimits(1, -0.3, 0.3);
      _circle_function_dn->SetParLimits(2, -0.1, 0.1);
      track_gr->Fit("circle_dn", "Q");
      double chi_dn = track_gr->GetFunction("circle_dn")->GetChisquare();
      _circle_function_up->SetParameters(80., 0., 0.);
      //_circle_function_up->SetParLimits(0, 5., 1000.);
      _circle_function_up->SetParLimits(1, -0.3, 0.3);
      _circle_function_up->SetParLimits(2, -0.1, 0.1);
      track_gr->Fit("circle_up", "Q");
      double chi_up = track_gr->GetFunction("circle_up")->GetChisquare();

      double best_chi, best_q;
      if (chi_up > chi_dn) {
        func = "circle_dn";
        best_chi = chi_dn;
        best_q = best_chi / track_gr->GetFunction("circle_up")->GetNDF();
      } else {
        func = "circle_up";
        best_chi = chi_up;
        best_q = best_chi / track_gr->GetFunction("circle_up")->GetNDF();
      }

      // TEMP fill comparison arc/lin fitting

      track_gr->Fit("pol1", "Q");
      lin_fit = track_gr->GetFunction("pol1");
      double lin_chi2 = lin_fit->GetChisquare();
      double lin_q = lin_chi2 / lin_fit->GetNDF();

      _chi2_ratio->Fill(best_chi / lin_chi2);
      _qulity_ratio->Fill(best_q / lin_q);
      // end of TEMP

      track_gr->Fit(func, "Q");
      fit = track_gr->GetFunction(func);
    }

    if (!fit)
      continue;

    float mom = fit->GetParameter(0) * units::B * units::clight / 1.e9;
    if (fit->GetParameter(2) < 0) mom *= -1.;
    _mom_reco->Fill(mom);
    _pos_reco->Fill(fit->GetParameter(2));
    _ang_reco->Fill(fit->GetParameter(1));

    double quality = fit->GetChisquare() / fit->GetNDF();

    if (_verbose > 2) {
      //fit = track_gr->GetFunction(func);

      std::cout << "Arc fit" << std::endl;
      std::cout << "radius\t" <<  fit->GetParameter(0) << std::endl;
      std::cout << "tan(a)\t" << fit->GetParameter(1) << std::endl;
      std::cout << "target\t" << 1/fit->GetParameter(2) << std::endl;
      std::cout << "q\t" << fit->GetChisquare() << "/" << fit->GetNDF() << std::endl;

      track_gr->Fit("pol1", "Q");
      fit = track_gr->GetFunction("pol1");
      std::cout << "Linear fit" << std::endl;
      std::cout << "k\t" << fit->GetParameter(1) << std::endl;
      std::cout << "b\t" << fit->GetParameter(0) << std::endl;
      std::cout << "q\t" << fit->GetChisquare() << "/" << fit->GetNDF() << std::endl;
    }

    _Chi2_track->Fill(quality);
    TF1* fit1[geom::nPadx];

    for (int i = 1; i < geom::nPadx-1; ++i) {
      if (!_do_arc_fit) {
        track_1[i]->Fit("pol1", "Q");
        fit1[i] = track_1[i]->GetFunction("pol1");
      } else {
        track_1[i]->Fit(func, "Q");
        fit1[i] = track_1[i]->GetFunction(func);
      }
      if (!fit1[i]) {
        std::cout << "Arc fit" << std::endl;
        std::cout << "r\t" <<  fit->GetParameter(0) << std::endl;
        std::cout << "x0\t" << fit->GetParameter(1) << std::endl;
        std::cout << "y0\t" << fit->GetParameter(2) << std::endl;
        std::cout << "q\t" << fit->GetChisquare() << "/" << fit->GetNDF() << std::endl;
      }
    }

    // second loop over columns
    for (auto row:track->GetCols()) {
      if (!row[0])
        continue;
      auto it_x = row[0]->GetCol();
      if (it_x == 0 || it_x == geom::nPadx-1)
        continue;

      if (track_pos[it_x]  == -999.)
        continue;

      double x    = geom::x_pos[it_x];
      double track_fit_y    = fit->Eval(x);
      double track_fit_y1   = fit1[it_x]->Eval(x);

      // fill SR
      _resol_col_hist[it_x]->Fill(track_pos[it_x] - track_fit_y);
      _resol_col_hist_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);

      if (cluster_N[it_x] == 2) {
        _resol_col_hist[it_x]->Fill(track_pos[it_x] - track_fit_y);
        _resol_col_hist_2pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      } else if (cluster_N[it_x] == 3) {
        _resol_col_hist_3pad[it_x]->Fill(track_pos[it_x] - track_fit_y);
        _resol_col_hist_3pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      }

      if (cluster_N[it_x] == 1)
        continue;
      // Fill PRF
      for (auto pad:row) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow();
        auto q = pad->GetQ();

        if (!cluster[it_x] || !q)
          continue;

        //double charge = 1. * q / cluster[it_x];
        double center_pad_y = geom::y_pos[it_y];

        // fill PRF
        _PRF_histo->Fill(center_pad_y - track_fit_y, q / a_peak[it_x]);
        _PRF_histo_col[it_x]->Fill(center_pad_y - track_fit_y, q / a_peak[it_x]);

        if (_verbose > 3)
          std::cout << "PRF fill " << center_pad_y - track_fit_y << "\t" << q / a_peak[it_x] << "\t( " << q << " / " <<  a_peak[it_x] << " )" << std::endl;

        if (cluster_N[it_x] == 2)
          _PRF_histo_2pad->Fill(center_pad_y - track_fit_y, q / a_peak[it_x]);
        else if (cluster_N[it_x] == 3)
          _PRF_histo_3pad->Fill(center_pad_y - track_fit_y, q / a_peak[it_x]);
        else if (cluster_N[it_x] == 4)
          _PRF_histo_4pad->Fill(center_pad_y - track_fit_y, q / a_peak[it_x]);
      }
    } // loop over colums
    delete track_gr;
    delete track_m;
    for (int i = 0; i < geom::nPadx; ++i) {
      delete track_1[i];
    }
  } // loop over tracks

  if (_store_event)
    _passed_events.push_back(event->GetID());

  return true;
}

bool SpatialResolAna::WriteOutput() {
  if (!_file_out)
    return true;

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

  TH1F* resol = new TH1F("resol", "", 1000, 0., 0.001);
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

    resol->Fill(sqrt(sigma * sigma_ex));

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

  std::cout << "Resol\t" << 1.e6*resol->GetMean() << " um" << "\tRMS\t" << 1.e6*resol->GetRMS() << " um" << std::endl;

  // Write objects
  AnalysisBase::WriteOutput();

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
  if (!ana->Initialize())               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}
