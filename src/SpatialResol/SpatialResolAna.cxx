#include "SpatialResolAna.hxx"

#include "Math/Functor.h"
#include "Fit/Fitter.h"

//! verbosity level
//! 1 (default) print progress, memory usage, time consumption
//! 2           print event and track number
//! 3           print fit details
//! 4           print PRF details

//********************************************************************
SpatialResolAna::SpatialResolAna(int argc, char** argv):
  AnalysisBase(argc, argv),
  _do_arc_fit(true),
  _do_full_track_fit(false),
  _do_separate_pad_fit(true),
  // WARNING
  _correction(false),
  _gaussian_residuals(true),
  _charge_uncertainty(true),
  _iteration(-1),
  _invert(false) {
  //********************************************************************

  // read CLI
  optind = 1;
  for (;;) {
    int c = getopt(argc, argv, "i:o:bv:drmst:ca");
    if (c < 0) break;
    switch (c) {
      case 't' : _iteration        = atoi(optarg);  break;
      case 'c' : _correction       = false;         break;
      case 'a' : _invert           = true;          break;
    }
  }

  if (_iteration == -1) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna(). Iteration should be defined as a input param" << std::endl;
    exit(1);
  }

  if (_do_separate_pad_fit && _do_full_track_fit) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna(). Incorrect input params definintion" << std::endl;
    exit(1);
  }

}

//********************************************************************
bool SpatialResolAna::Initialize() {
  //********************************************************************
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
    prev_file_name   = prev_file_name(0, prev_file_name.Index("iter"));
    prev_file_name  += "iter";
    prev_file_name  += TString::Itoa(_iteration - 1, 10);
    prev_file_name  += ".root";

    _Prev_iter_file = new TFile(prev_file_name.Data(), "READ");
    _PRF_function   = (TF1*)_Prev_iter_file->Get("PRF_function");
    auto uncertainty_graph = (TGraphErrors*)_Prev_iter_file->Get("resol_final");
    _uncertainty = 0;
    for (auto i = 0; i < uncertainty_graph->GetN(); ++i) {
      double x, y;
      uncertainty_graph->GetPoint(i, x, y);
      _uncertainty += y / uncertainty_graph->GetN();
    }

    // read event list passed reconstruction+selection
    // at the moment skip for the 1 iteration (after 0) files with TEvent
    // reason: in at the step 0 TEvent file is generated need to look at all events in it
    if (!_work_with_event_file || _iteration != 1) {
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
/*
    _PRF_function = new TF1("PRF_function",
      " [0] * exp(-4*TMath::Log(2)*(1-[1])*TMath::Power(x/[2], 2.)) / (1+4 * [1] * TMath::Power(x/[2], 2.) )", prf_min, prf_max);
    _PRF_function->SetParName(0, "Const");
    _PRF_function->SetParName(1, "r");
    _PRF_function->SetParName(2, "w");

    auto c = 1.;
    auto r = 0.5;
    auto s = 0.005;
    _PRF_function->SetParameters(c, r, s);

*/
    _PRF_function  = new TF1("PRF_function", "[0]*(1+[1]*x*x + [2] * x*x*x*x) / (1+[3]*x*x+[4]*x*x*x*x)", prf_min, prf_max);
    _PRF_function->SetParName(0, "Const");
    _PRF_function->SetParName(1, "a2");
    _PRF_function->SetParName(2, "a4");
    _PRF_function->SetParName(3, "b2");
    _PRF_function->SetParName(4, "b4");

    double co = 1.;
    double a2 = 2.35167e3;
    double a4 = 6.78962e7;
    double b2 = 3.36748e3;
    double b4 = 6.45311e8;
    _PRF_function->SetParameters(co, a2, a4, b2, b4);

    if (_do_full_track_fit)
      _PRF_function->FixParameter(0, 1.);

  }

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

  for (auto i = 0; i < GetMaxColumn(); ++i)
    _PRF_histo_col[i] = new TH2F(Form("PRF_histo_col_%i", i),"", prf_bin, prf_min, prf_max, 150,0.,1.5);

  _PRF_graph = new TGraphErrors();
  _PRF_graph->SetName("PRF_graph");

  for (auto j = 0; j < GetMaxColumn(); ++j) {
    _resol_col_hist[j]  = new TH1F(Form("resol_histo_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_except[j]  = new TH1F(Form("resol_histo1__%i", j), "", resol_bin, resol_min, resol_max);

    _resol_col_hist_2pad[j]  = new TH1F(Form("resol_histo_2pad_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_2pad_except[j]  = new TH1F(Form("resol_histo1_2pad_%i", j), "", resol_bin, resol_min, resol_max);

    _resol_col_hist_3pad[j]  = new TH1F(Form("resol_histo_3pad_%i", j), "", resol_bin, resol_min, resol_max);
    _resol_col_hist_3pad_except[j]  = new TH1F(Form("resol_histo1_3pad_%i", j), "", resol_bin, resol_min, resol_max);
  }

  _residual_sigma_unbiased  = new TGraphErrors();
  _residual_sigma_biased    = new TGraphErrors();
  _residual_sigma           = new TGraphErrors();

  _residual_mean            = new TGraphErrors();

  _residual_sigma_unbiased->SetName("resol_unb");
  _residual_sigma_biased->SetName("resol_bia");
  _residual_sigma->SetName("resol_final");
  _residual_mean->SetName("mean");

  _Chi2_track = new TH1F("Chi2_Track", "", 1000, 0., 100.);
  _Cols_used  = new TH1F("cols_used", "", 40, 0., 40.);

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
  _output_vector.push_back(_Cols_used);
  _output_vector.push_back(_qulity_ratio);
  _output_vector.push_back(_chi2_ratio);
  _output_vector.push_back(_mom_reco);
  _output_vector.push_back(_pos_reco);
  _output_vector.push_back(_ang_reco);

  for (auto j = 0; j < GetMaxColumn(); ++j) {
    _output_vector.push_back(_resol_col_hist[j]);
    _output_vector.push_back(_resol_col_hist_except[j]);
  }

  _uncertainty_vs_prf_gr = new TGraphErrors();
  _output_vector.push_back(_uncertainty_vs_prf_gr);

  if (_do_separate_pad_fit) {
    Double_t arr[4] = {0., .15, .7, 1.01};
    _prf_scale_axis = new TAxis(3, arr);
    for (auto j = 0; j < GetMaxColumn(); ++j) {
      for (auto id=0; id < 3; ++id) {
        _Fit_quality_plots[id][j] = new TH1F(Form("pad_fit_q_%i_%i", j, id), "", 300, 0., 30.);
        _output_vector.push_back(_Fit_quality_plots[id][j]);
      }
    }
  }

  _x_scan_axis = new TAxis(x_scan_bin, x_scan_min, x_scan_max);
  for (auto j = 0; j < GetMaxColumn(); ++j) {
    for (auto i = 0; i < x_scan_bin; ++i) {
      _resol_col_x_scan[j][i] = new TH1F(Form("resol_histo_Xscan_%i_%i", j, i), "", resol_bin, resol_min, resol_max);
      _output_vector.push_back(_resol_col_x_scan[j][i]);
      _mult_x_scan[j][i] = new TH1F(Form("mult_histo_Xscan_%i_%i", j, i), "multiplicity", 10, 0., 10.);
      _output_vector.push_back(_mult_x_scan[j][i]);
    }
  }

  for (auto i = 0; i < prf_error_bins; ++i) {
    _uncertainty_prf_bins[i] = new TH1F(Form("error_prf_bin_%i", i), "", resol_bin, resol_min, resol_max);
    _output_vector.push_back(_uncertainty_prf_bins[i]);
  }

  for (auto j = 0; j < GetMaxColumn(); ++j) {
    for (auto i = 0; i < x_scan_bin; ++i) {
      _resol_col_x_scan_lim_mult[j][i] = new TH1F(Form("resol_histo_Xscan_%i_%i", j, i), "", resol_bin, resol_min, resol_max);
      _output_vector.push_back(_resol_col_x_scan_lim_mult[j][i]);
    }
  }

  for (auto j = 0; j < GetMaxColumn(); ++j) {
    _output_vector.push_back(_PRF_histo_col[j]);
  }

  _passed_events.clear();

  std::cout << "done" << std::endl;
  std::cout << "      PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula();
  std::cout << "  with ";
  for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
    std::cout << "  " << _PRF_function->GetParameter(i) << ",";
  std::cout << std::endl;

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  // Initialize timers
  _sw_partial[2] = new TStopwatch();
  _sw_partial[2]->Reset();
  _sw_partial[3] = new TStopwatch();
  _sw_partial[3]->Reset();
  _sw_partial[4] = new TStopwatch();
  _sw_partial[4]->Reset();

  return true;
}

//********************************************************************
double SpatialResolAna::GetTrackPosInPad(const std::vector<THit*>& col,
    const int cluster, const double pos,
    std::vector<std::vector<std::pair< double, std::pair<double, double> > > >& pos_in_pad) {
  //********************************************************************
  if (!_iteration) {
    pos_in_pad[col[0]->GetCol(_invert)].push_back(std::make_pair(0.1, std::make_pair(pos, default_error)));
    return pos;
  }

  double sum1 = 0;
  double sum2 = 0;

  double track_pos;
  double track_pos_err;

  for (auto pad:col) {
    auto q      = pad->GetQ();
    auto it_y   = pad->GetRow(_invert);
    auto it_x   = pad->GetCol(_invert);
    if (!q)
      continue;

    if (_verbose > 5)
      std::cout << "pad " << "\t" << 1.*q/cluster << "\t" << _PRF_function->Eval(prf_min) << std::endl;

    double center_pad_y = GetYpos(it_y, _invert);
    if (1.*q/cluster < _PRF_function->Eval(fit_bound_left))
      continue;

    if (_verbose > 5)
      std::cout << q << "\t" << cluster << std::endl;

    if (1.*q/cluster > _PRF_function->GetParameter(0)) {
      track_pos = center_pad_y;
      track_pos_err = 0.003;
    } else {
      double pos_bias = _PRF_function->GetX(1.*q/cluster, 0., fit_bound_right);

      if (pos > center_pad_y)
        track_pos = center_pad_y + pos_bias;
      else
        track_pos = center_pad_y - pos_bias;

      track_pos_err = sigma_pedestal / cluster;
      track_pos_err /= abs(_PRF_function->Derivative(pos_bias));

      if (_verbose > 5)
        std::cout << "errors\t" << sigma_pedestal / cluster << "\t" << _PRF_function->Derivative(pos_bias) << "\t" << track_pos_err << std::endl;
    }

    if (1.*q/cluster > 0.5 && track_pos_err > 0.003)
      track_pos_err = 0.003;

    if (track_pos_err > 0.08)
      track_pos_err = 0.08;

    if (_verbose > 5)
      std::cout << "pad pos \t" << track_pos << "\t" << track_pos_err << std::endl;

    pos_in_pad[it_x].push_back(std::make_pair(1.*q/cluster, std::make_pair(track_pos, track_pos_err)));

    sum1 += track_pos * TMath::Power(track_pos_err, -2);
    sum2 += TMath::Power(track_pos_err, -2);
  } // loop over pads in cluster

  return sum1 / sum2;
}


//********************************************************************
double SpatialResolAna::GetClusterPosCERN(const std::vector<THit*>& col,
    const int cluster, const double pos) {
  //********************************************************************

  if (!_iteration)
    return pos;

  auto chi2Function_cluster = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;

      double a = 1. * q / cluster;
      double center_pad_y = GetYpos(it_y, _invert);

      // avoid using pads wich are far away from track
      // limit by PRF fitting range (PRF function robustness)
      if (abs(center_pad_y - pos) > fit_bound_right)
        continue;

      double part = (a - _PRF_function->Eval(par[0] - center_pad_y));
      if (_charge_uncertainty) {
        double c = 1.*cluster;
        double b = 1.*q;
        part *= c*c;
        part /= c*sqrt(b) + b*sqrt(c);
      }
      part *= part;

      chi2 += part;
    }
    return chi2;
  };

  ROOT::Math::Functor fcn_cluster(chi2Function_cluster,1);
  ROOT::Fit::Fitter  fitter_cluster;

  double pStart[1] = {pos};
  fitter_cluster.SetFCN(fcn_cluster, pStart);
  fitter_cluster.Config().ParSettings(0).SetName("y");

  bool ok = fitter_cluster.FitFCN();
  (void)ok;
  const ROOT::Fit::FitResult & result_cluster = fitter_cluster.Result();

  return result_cluster.GetParams()[0];
}

// BUG PRF limitation is not implemented in the ILC fitting procedure
//********************************************************************
double SpatialResolAna::GetClusterPosILC(const std::vector<THit*>& col,
  const double pos) {
  //********************************************************************

  auto chi2Function_cluster = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;

    double a_nom = 0.;
    double a_den = 0.;
    double a_tot = 0.;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;
      double center_pad_y = GetYpos(it_y, _invert);

      a_nom += _PRF_function->Eval(center_pad_y - par[0]);
      a_den += TMath::Power(_PRF_function->Eval(center_pad_y - par[0]), 2) / q;
    }
    a_tot = a_nom / a_den;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;
      double center_pad_y = GetYpos(it_y, _invert);

      double part = (q - a_tot*_PRF_function->Eval(center_pad_y - par[0]));
      part *= part;
      part /= q;

      chi2 += part;
    }

    return chi2;
  };

  ROOT::Math::Functor fcn_cluster(chi2Function_cluster,1);
  ROOT::Fit::Fitter  fitter_cluster;

  double pStart[1] = {pos};
  fitter_cluster.SetFCN(fcn_cluster, pStart);
  fitter_cluster.Config().ParSettings(0).SetName("y");

  bool ok = fitter_cluster.FitFCN();
  (void)ok;
  const ROOT::Fit::FitResult & result_cluster = fitter_cluster.Result();

  return result_cluster.GetParams()[0];
}

//********************************************************************
TF1* SpatialResolAna::GetTrackFitCERN(const double* track_pos,
  const int* mult, const int miss_id) {
  //********************************************************************
  TGraphErrors* track_gr = new TGraphErrors();

  for (auto it_x = 0; it_x < GetMaxColumn(); ++it_x) {
    double x = GetXpos(it_x, _invert);
    if (track_pos[it_x] == -999.)
      continue;

    if (it_x == miss_id)
      continue;

    track_gr->SetPoint(track_gr->GetN(), x, track_pos[it_x]);
    double error;
    if (mult[it_x] == 1)
      error = one_pad_error;
    else {
      if (_iteration == 0)
        error = default_error;
      else
        error = _uncertainty;
    }
    track_gr->SetPointError(track_gr->GetN()-1, 0., error);
  } // loop over x

  TF1* fit;
  TString func;
  if (!_do_arc_fit) {
    track_gr->Fit("pol1", "Q");
    fit = (TF1*)track_gr->GetFunction("pol1")->Clone();
    func = "pol1";
  } else {
    Float_t q_up, q_down;
    q_up = q_down = 1.e9;
    _circle_function_dn->SetParameters(80., 0, 0.);
    track_gr->Fit("circle_dn", "Q");
    fit = track_gr->GetFunction("circle_dn");
    if (fit)
      q_down = fit->GetChisquare() / fit->GetNDF();

    _circle_function_up->SetParameters(80., 0, 0.);
    track_gr->Fit("circle_up", "Q");
    fit = track_gr->GetFunction("circle_up");
    if (fit)
      q_up = fit->GetChisquare() / fit->GetNDF();

    if (q_up > q_down)
      func = "circle_dn";
    else
      func = "circle_up";

    track_gr->Fit(func, "Q");
    if (!track_gr->GetFunction(func))
      return NULL;
    fit = (TF1*)track_gr->GetFunction(func)->Clone();

  }

  delete track_gr;

  if (!fit)
    return NULL;

  double quality = fit->GetChisquare() / fit->GetNDF();
  _Chi2_track->Fill(quality);

  return fit;
}

//********************************************************************
TF1* SpatialResolAna::GetTrackFitSeparatePad(const std::vector<std::vector<std::pair< double, std::pair<double, double> > > >pos_in_pad,
  const int miss_id) {
  //********************************************************************
  TGraphErrors* track_gr = new TGraphErrors();

  for (auto it_x = 0; it_x < GetMaxColumn(); ++it_x) {
    double x = GetXpos(it_x, _invert);
    if (!pos_in_pad[it_x].size())
      continue;

    if (it_x == miss_id)
      continue;

    for (auto it = pos_in_pad[it_x].begin(); it < pos_in_pad[it_x].end(); ++it) {
      track_gr->SetPoint(track_gr->GetN(), x, (*it).second.first);
      track_gr->SetPointError(track_gr->GetN()-1, 0., (*it).second.second);
    }
  } // loop over x

  TF1* fit;
  TString func;
  if (!_do_arc_fit) {
    track_gr->Fit("pol1", "Q");
    fit = (TF1*)track_gr->GetFunction("pol1")->Clone();
    func = "pol1";
  } else {
    Float_t q_up, q_down;
    q_up = q_down = 1.e9;
    _circle_function_dn->SetParameters(80., 0, 0.);
    track_gr->Fit("circle_dn", "Q");
    fit = track_gr->GetFunction("circle_dn");
    if (fit)
      q_down = fit->GetChisquare() / fit->GetNDF();

    _circle_function_up->SetParameters(80., 0, 0.);
    track_gr->Fit("circle_up", "Q");
    fit = track_gr->GetFunction("circle_up");
    if (fit)
      q_up = fit->GetChisquare() / fit->GetNDF();

    if (q_up > q_down)
      func = "circle_dn";
    else
      func = "circle_up";

    track_gr->Fit(func, "Q");
    if (!track_gr->GetFunction(func))
      return NULL;
    fit = (TF1*)track_gr->GetFunction(func)->Clone();

  }

  delete track_gr;

  if (!fit)
    return NULL;

  double quality = fit->GetChisquare() / fit->GetNDF();
  _Chi2_track->Fill(quality);

  return fit;
}

//********************************************************************
bool SpatialResolAna::MissColumn(const int it_x) {
  //********************************************************************

  if (it_x == 0 || it_x == GetMaxColumn()-1)
    return true;

  return false;
}

Int_t SpatialResolAna::GetMaxColumn() {
  if (!_invert)
    return geom::nPadx;
  else
    return geom::nPady;
}

// TODO move it to selection
//********************************************************************
bool SpatialResolAna::UseCluster(const std::vector<THit*>& col) {
  //********************************************************************
  return true;
  /*auto N = std::count_if(col.begin(), col.end(),
                       [](const THit* h1){return (h1->GetQ() > 0); });

  (void)N;

  //if (N != 3)
  //  return false;

  //return true;
  auto it_max = std::max_element(col.begin(), col.end(),
                       [](const THit* h1, const THit* h2) { return h1->GetQ() < h2->GetQ(); });
  int maxQ = (it_max == col.end()) ? 0 : (*it_max)->GetQ();

  if (maxQ > 2000)
    return false;

  return true;*/
}

//********************************************************************
TF1* SpatialResolAna::GetTrackFitILC(const TTrack* track, const double pos,
  const int miss_id) {
  //********************************************************************
  auto chi2Function = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;
    float y_pos = 0;

    if (par[3] > 0)
      _circle_function_up->SetParameters(par[0], par[1], par[2]);
    else
      _circle_function_dn->SetParameters(par[0], par[1], par[2]);


    for (auto col:track->GetCols(_invert)) {
      auto it_x = col[0]->GetCol(_invert);
      auto x    = GetXpos(it_x, _invert);

      // ommit first and last
      if (MissColumn(it_x))
        continue;

      if (!UseCluster(col))
        continue;

      if (it_x == miss_id)
        continue;

      double a_nom = 0.;
      double a_den = 0.;
      double a_tot = 0.;

      if (par[3] > 0) {
        y_pos = _circle_function_up->Eval(x);
      } else {
        y_pos = _circle_function_dn->Eval(x);
      }

      for (auto pad:col) {
        auto q      = pad->GetQ();
        auto it_y   = pad->GetRow(_invert);
        if (!q)
          continue;

        double center_pad_y = GetYpos(it_y, _invert);

        a_nom += _PRF_function->Eval(center_pad_y - y_pos);
        a_den += TMath::Power(_PRF_function->Eval(center_pad_y - y_pos), 2) / q;
      }
      a_tot = a_nom / a_den;

      for (auto pad:col) {
        auto q      = pad->GetQ();
        auto it_y   = pad->GetRow(_invert);
        if (!q)
          continue;
        double center_pad_y = GetYpos(it_y, _invert);

        double part = (q - a_tot*_PRF_function->Eval(center_pad_y - y_pos));
        part *= part;
        part /= q;

        chi2 += part;
      }
    }
    return chi2;
  };

  float quality_dn, quality_up;
  quality_dn = quality_up = 1.e9;

  _sw_partial[2]->Stop();
  _sw_partial[3]->Start(false);

  ROOT::Math::Functor fcn(chi2Function,4);
  ROOT::Fit::Fitter  fitter_dn;

  // fitting arc down
  double pStart[4] = {80., 0., pos, -1};
  fitter_dn.SetFCN(fcn, pStart, 4, true);
  fitter_dn.Config().ParSettings(0).SetName("R");
  fitter_dn.Config().ParSettings(1).SetName("sin(#alpha)");
  fitter_dn.Config().ParSettings(2).SetName("Target");
  fitter_dn.Config().ParSettings(3).SetName("Arc_dir (up/down)");

  fitter_dn.Config().ParSettings(0).SetLimits(0., 1.e5);
  fitter_dn.Config().ParSettings(3).Fix();

  bool ok = fitter_dn.FitFCN();
  (void)ok;
  const ROOT::Fit::FitResult & result_dn = fitter_dn.Result();

  _circle_function_dn->SetParameters( result_dn.GetParams()[0],
                                      result_dn.GetParams()[1],
                                      result_dn.GetParams()[2]);
  if (_verbose > 2)
    result_dn.Print(std::cout);

  quality_dn = result_dn.Chi2() / (track->GetHits().size() - 3);

  // fitting arc up
  ROOT::Math::Functor fcn_up(chi2Function,4);
  ROOT::Fit::Fitter  fitter_up;
  double pStart_up[4] = {80., 0., pos, 1};
  fitter_up.SetFCN(fcn_up, pStart_up, 4, true);
  fitter_up.Config().ParSettings(0).SetName("R");
  fitter_up.Config().ParSettings(1).SetName("sin(#alpha)");
  fitter_up.Config().ParSettings(2).SetName("Target");
  fitter_up.Config().ParSettings(3).SetName("Arc_dir (up/down)");
  fitter_up.Config().ParSettings(0).SetLimits(0., 1.e5);
  fitter_up.Config().ParSettings(3).Fix();

  ok = fitter_up.FitFCN();

  const ROOT::Fit::FitResult & result_up = fitter_up.Result();

  _circle_function_up->SetParameters( result_up.GetParams()[0],
                                      result_up.GetParams()[1],
                                      result_up.GetParams()[2]);
  if (_verbose > 2)
    result_up.Print(std::cout);

  quality_up = result_up.Chi2() / (track->GetHits().size() - 3);

  TF1* fit;
  TString func;
  float quality;
  if (quality_up > quality_dn) {
    fit = _circle_function_dn;
    quality = quality_dn;
  } else {
    fit = _circle_function_up;
    quality = quality_up;
  }

  if (quality > 100.)
    return NULL;

  _Chi2_track->Fill(quality);

  return fit;
}

//********************************************************************
bool SpatialResolAna::ProcessEvent(const TEvent* event) {
  //********************************************************************

  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    if (_verbose > 1)
      std::cout << "Track id = " << trackId << std::endl;

    if(sel::GetNonZeroCols(track, _invert).size() != (uint)GetMaxColumn()) return false;
    if(sel::GetColsMaxSep(track, _invert)>8) return false;
    if (sel::GetColsMaxGap(track, _invert) > 0) return false;
    std::vector<double> fit_v = sel::GetFitParams(track, _invert);
    if(fit_v[0]>1.0e6) return false;
    if (abs(fit_v[2]) > 0.08) return false;

    _store_event = true;

    if (_verbose > 1)
      std::cout << "selected" << std::endl;

    int cluster[geom::nPadx];
    int cluster_N[geom::nPadx];
    double track_pos[geom::nPadx];
    double cluster_mean[geom::nPadx];
    float charge_max[geom::nPadx];
    double a_peak_fit[geom::nPadx];

    std::vector<std::vector<std::pair<double, std::pair<double, double> > > >pos_in_pad;
    pos_in_pad.clear();
    pos_in_pad.resize(GetMaxColumn());
    int Ndots = 0;

    // At the moment ommit first and last column
    // first loop over columns
    _sw_partial[2]->Start(false);
    for (auto col:track->GetCols(_invert)) {
      if (!col[0])
        continue;
      auto it_x = col[0]->GetCol(_invert);

      cluster[it_x]       = 0;
      cluster_N[it_x]     = 0;
      charge_max[it_x]    = 0;
      track_pos[it_x]     = -999.;
      cluster_mean[it_x]  = 0.;
      a_peak_fit[it_x]    = 0.;

      // exlude 1st/last column
      if (MissColumn(it_x))
        continue;

      TH1F* cluster_h;
      if (!_invert)
        cluster_h = new TH1F("cluster", "", geom::nPady, -1.*geom::MM_dy - geom::dy/2., geom::MM_dy + geom::dy/2.);
      else
        cluster_h = new TH1F("cluster", "", geom::nPadx, -1.*geom::MM_dx - geom::dx/2., geom::MM_dx + geom::dx/2.);

      // loop over rows
      for (auto pad:col) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow(_invert);
        auto q = pad->GetQ();

        cluster[it_x] += q;
        ++cluster_N[it_x];
        cluster_h->Fill(GetYpos(it_y, _invert), q);
        if (charge_max[it_x] < q)
          charge_max[it_x] = q;
      } // end of loop over rows

      cluster_mean[it_x] = cluster_h->GetMean();

      delete cluster_h;
      cluster_h = NULL;

      if (!cluster[it_x] || !charge_max[it_x])
        continue;

      // cluster usage cut. e.g. max charge
      if (!UseCluster(col))
        continue;

      ++Ndots;

      if (_do_full_track_fit)
        track_pos[it_x] = GetClusterPosILC(col, cluster_mean[it_x]);
      else if (_do_separate_pad_fit)
        track_pos[it_x] = GetTrackPosInPad(col,
          cluster[it_x], cluster_mean[it_x], pos_in_pad);
      else
        track_pos[it_x] = GetClusterPosCERN(col,
          cluster[it_x], cluster_mean[it_x]);

      if (_verbose > 5)
        std::cout << "Cluster pos " << track_pos[it_x] << "\t" << cluster_mean[it_x] << std::endl;
    } // loop over columns

    if (_verbose > 3)
      std::cout << "Loop over columns done" << std::endl;

    _Cols_used->Fill(Ndots);

    _sw_partial[2]->Stop();
    _sw_partial[3]->Start(false);

    TF1* fit;
    if (_do_full_track_fit) {
      fit = GetTrackFitILC(track, track_pos[1]);
      if (fit)
        fit = (TF1*)GetTrackFitILC(track, track_pos[1])->Clone();
      //fit = (TF1*)GetTrackFitILC(track, track_pos[1])->Clone();
    } else if (_do_separate_pad_fit) {
      fit = GetTrackFitSeparatePad(pos_in_pad);
    } else
      fit = GetTrackFitCERN(track_pos, cluster_N);

    if (_verbose > 3)
      std::cout << "Track fit done" << std::endl;

    if (!fit)
      continue;

    TString func = fit->GetName();

    // in case of arc fitting fill the momentum histo
    if (_do_arc_fit) {
      float mom = fit->GetParameter(0) * units::B * units::clight / 1.e9;
      if (func.CompareTo("circle_dn") == 0) mom *= -1.;
      _mom_reco->Fill(mom);
      _pos_reco->Fill(fit->GetParameter(2));
      _ang_reco->Fill(fit->GetParameter(1));
    }

    TF1* fit1[geom::nPadx];
    for (int i = 1; i < GetMaxColumn()-1; ++i) {
      if (!_correction)
        fit1[i] = fit;
      else {
        if (_do_full_track_fit)
          fit1[i] = (TF1*)GetTrackFitILC(track, track_pos[1], i)->Clone();
        else if (_do_separate_pad_fit)
          fit1[i] = GetTrackFitSeparatePad(pos_in_pad, i);
        else
          fit1[i] = GetTrackFitCERN(track_pos, cluster_N, i);

      }
    }

    _sw_partial[3]->Stop();
    _sw_partial[4]->Start(false);

    // second loop over columns
    for (auto col:track->GetCols(_invert)) {
      if (!col[0])
        continue;
      auto it_x = col[0]->GetCol(_invert);
      if (MissColumn(it_x))
        continue;

      if (track_pos[it_x]  == -999.)
        continue;

      double x    = GetXpos(it_x, _invert);
      double track_fit_y    = fit->Eval(x);
      double track_fit_y1   = fit1[it_x]->Eval(x);

      if (_verbose > 4) {
        std::cout << "Residuals cluster:track\t" << track_pos[it_x] << "\t" << track_fit_y << "\t" << (track_pos[it_x] - track_fit_y)*1e6 << std::endl;
      }

      // fill SR
      _resol_col_hist[it_x]->Fill(track_pos[it_x] - track_fit_y);
      _resol_col_hist_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);

      // fill x scan
      auto bin = _x_scan_axis->FindBin(track_fit_y);

      if (bin > 0 && bin <= x_scan_bin) {
        _resol_col_x_scan[it_x][bin-1]->Fill(track_pos[it_x] - track_fit_y);
        _mult_x_scan[it_x][bin-1]->Fill(cluster_N[it_x]);

        if (cluster_N[it_x] == 3 || cluster_N[it_x] == 4)
          _resol_col_x_scan_lim_mult[it_x][bin-1]->Fill(track_pos[it_x] - track_fit_y);
      }

      if (cluster_N[it_x] == 2) {
        _resol_col_hist[it_x]->Fill(track_pos[it_x] - track_fit_y);
        _resol_col_hist_2pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      } else if (cluster_N[it_x] == 3) {
        _resol_col_hist_3pad[it_x]->Fill(track_pos[it_x] - track_fit_y);
        _resol_col_hist_3pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      }

      // in case of full track fit calc normalisation coefficient
      if (_do_full_track_fit) {
        float a_nom = 0.;
        float a_den = 0.;
        for (auto pad:col) {
          auto q      = pad->GetQ();
          auto it_y   = pad->GetRow(_invert);
          if (!q)
            continue;
          double center_pad_y = GetYpos(it_y, _invert);

          if (abs(center_pad_y - track_fit_y) > fit_bound_right)
            continue;

          a_nom += _PRF_function->Eval(center_pad_y - track_fit_y);
          a_den += TMath::Power(_PRF_function->Eval(center_pad_y - track_fit_y), 2) /   q;
        }

        a_peak_fit[it_x] = a_nom / a_den;
      } else
        a_peak_fit[it_x] = 1.*cluster[it_x];

      // fill pad accurace
      if (_do_separate_pad_fit) {
        for (auto it = pos_in_pad[it_x].begin();
                  it < pos_in_pad[it_x].end();
                  ++it) {
          int bin = -1 + _prf_scale_axis->FindBin((*it).first);
        // TODO make the definition through constants
          if (bin < 0 || bin > 2)
            std::cout << "Error bin " << bin << "\t" << (*it).first << std::endl;
          _Fit_quality_plots[bin][it_x]->Fill(abs(track_fit_y - (*it).second.first) / (*it).second.second);

          bin = -1 + _prf_error_axis->FindBin((*it).first);
          if (bin < 0 || bin > prf_error_bins)
            std::cout << "Error bin " << bin << "\t" << (*it).first << std::endl;

          _uncertainty_prf_bins[bin]->Fill(track_fit_y - (*it).second.first);
        }
      }

      if (cluster_N[it_x] == 1)
        continue;
      // Fill PRF
      for (auto pad:col) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow(_invert);
        auto q = pad->GetQ();

        if (!cluster[it_x] || !q)
          continue;

        double center_pad_y = GetYpos(it_y, _invert);

        // fill PRF
        _PRF_histo->Fill( center_pad_y - track_fit_y,
                          q / a_peak_fit[it_x]);
        _PRF_histo_col[it_x]->Fill( center_pad_y - track_fit_y,
                                    q / a_peak_fit[it_x]);

        if (cluster_N[it_x] == 2)
          _PRF_histo_2pad->Fill(center_pad_y - track_fit_y,
                                q / a_peak_fit[it_x]);
        else if (cluster_N[it_x] == 3)
          _PRF_histo_3pad->Fill(center_pad_y - track_fit_y,
                                q / a_peak_fit[it_x]);
        else if (cluster_N[it_x] == 4)
          _PRF_histo_4pad->Fill(center_pad_y - track_fit_y,
                                q / a_peak_fit[it_x]);
      }
    } // loop over colums
    _sw_partial[4]->Stop();
    for (int i = 0; i < GetMaxColumn(); ++i)
      if (fit1[i] && _correction) {
        delete fit1[i];
        fit1[i] = NULL;
      }
    //  delete track_1[i];
    //}
    delete fit;

    if(_test_mode) this->DrawSelectionCan(event,trackId);
  } // loop over tracks

  if (_store_event)
    _passed_events.push_back(event->GetID());

  return true;
}

TCanvas* SpatialResolAna::DrawSelectionCan(const TEvent* event, int trkID) {
  TH2F    *MM      = new TH2F("MM","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TH2F    *MMsel   = new TH2F("MMsel","", geom::nPadx,geom::x_pos[0] - geom::dx, geom::x_pos[geom::nPadx-1]+geom::dx,
                                          geom::nPady,geom::y_pos[0] - geom::dy, geom::y_pos[geom::nPady-1]+geom::dy);
  TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");

  for (auto x = 0; x < geom::nPadx; ++x) {
    for (auto y = 0; y < geom::nPady; ++y) {
      auto max = 0;
      for (auto t = 0; t < geom::Nsamples; ++t) {
        if (_padAmpl[x][y][t] > max) {
          max = _padAmpl[x][y][t];
        }
      } // over t
      if (max)
        MM->Fill(x, y, max);
    }
  }

  // sel hits
  for (auto h:event->GetTracks()[trkID]->GetHits()){
    if(!h->GetQ()) continue;
    event3D->Fill(h->GetTime(),h->GetRow(),h->GetCol(),h->GetQ());
    MMsel->Fill(geom::y_pos[h->GetCol()],geom::x_pos[h->GetRow()],h->GetQ());
  }

  TCanvas *canv = new TCanvas("canv", "canv", 0., 0., 1400., 600.);
  canv->Divide(3,1);
  canv->cd(1);
  if (MM->Integral())
    MM->Draw("colz");
  //_PRF_histo->Draw("COLZ");
  canv->cd(2);
  MMsel->Draw("COLZ");

  canv->cd(3);
  event3D->Draw("x:y:z:c","","box2");
  TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetLimits(0,geom::nPadx);
  htemp->GetYaxis()->SetLimits(0,geom::nPady);
  htemp->GetZaxis()->SetLimits(0,500);
  htemp->SetTitle("");
  canv->Update();
  return canv;
  //canv->WaitPrimitive();
  /*delete htemp;
  delete canv;

  delete MM;
  delete MMsel;
  delete event3D;*/
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
  for (auto i = 1; i < GetMaxColumn() - 1; ++i) {

    TH1F* res     = _resol_col_hist[i];
    TH1F* res_e   =  _resol_col_hist_except[i];

    Float_t mean, sigma, sigma_ex;
    Float_t mean_e, sigma_e, sigma_ex_e;
    mean_e = sigma_e = sigma_ex_e = 0.;

    mean = res->GetMean();
    mean_e = res->GetMeanError();
    float max   = res->GetMaximum();
    float start = res->GetBinLowEdge(res->FindFirstBinAbove(max/2));
    float end   = res->GetBinLowEdge(res->FindLastBinAbove(max/2)) +
    res->GetBinWidth(res->FindLastBinAbove(max/2));

    sigma = 0.5 * (end - start);

    if (!res->Integral()) {
      std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Empty residuals" << std::endl;
      exit(1);
    }

    if (_gaussian_residuals) {
      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;
      _resol_col_hist_2pad[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad[i]->Fit("gaus", "Q");

      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        exit(1);
      }

      res_e->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;
      _resol_col_hist_2pad_except[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad_except[i]->Fit("gaus", "Q");

      mean     = func->GetParameter(1);
      sigma    = func->GetParameter(2);
      sigma_ex =res_e->GetFunction("gaus")->GetParameter(2);

      mean_e      = func->GetParError(1);
      sigma_e     = func->GetParError(2);
      sigma_ex_e  = res_e->GetFunction("gaus")->GetParError(2);
    } else {
      // use FWHM
      max   = res_e->GetMaximum();
      start = res_e->GetBinLowEdge(res_e->FindFirstBinAbove(max/2));
      end   = res_e->GetBinLowEdge(res_e->FindLastBinAbove(max/2)) +
      res_e->GetBinWidth(res_e->FindLastBinAbove(max/2));

      sigma_ex = 0.5 * (end - start);
    }

    _residual_sigma_biased->SetPoint(_residual_sigma_biased->GetN(),
                                      i+1, sigma);
    _residual_sigma_biased->SetPointError(_residual_sigma_biased->GetN()-1,
                                      0, sigma_e);

    _residual_sigma_unbiased->SetPoint(_residual_sigma_unbiased->GetN(),
                                      i+1, sigma_ex);
    _residual_sigma_unbiased->SetPointError(_residual_sigma_unbiased->GetN()-1,
                                      0, sigma_ex);

    Float_t val, err;
    val = sqrt(sigma * sigma_ex);
    if (val)
      err = 0.5 / val;
    else err = 0.;
    err *= sigma_e * sigma_ex + sigma * sigma_ex_e;
    _residual_sigma->SetPoint(_residual_sigma->GetN(), i+1, val);
    _residual_sigma->SetPointError(_residual_sigma->GetN()-1, 0, err);

    _residual_mean->SetPoint(_residual_mean->GetN(), i+1, mean);
    _residual_mean->SetPointError(_residual_mean->GetN()-1, 0, mean_e);

    // print out monitoring function
    resol->Fill(sqrt(sigma * sigma_ex));
  } // loop over column

  _PRF_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);

  if (_do_separate_pad_fit) {
    for (auto prf_bin = 0; prf_bin < prf_error_bins; ++prf_bin) {
      TH1F* res = _uncertainty_prf_bins[prf_bin];
      Float_t mean, sigma, sigma_e;
      mean = res->GetMean();
      float max   = res->GetMaximum();
      float start = res->GetBinLowEdge(res->FindFirstBinAbove(max/2));
      float end   = res->GetBinLowEdge(res->FindLastBinAbove(max/2)) +
      res->GetBinWidth(res->FindLastBinAbove(max/2));

      sigma = 0.5 * (end - start);

      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;

      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        exit(1);
      }

      mean     = func->GetParameter(1);
      sigma    = func->GetParameter(2);
      sigma_e  = func->GetParError(2);

      _uncertainty_vs_prf_gr->SetPoint(_uncertainty_vs_prf_gr->GetN(),
                              _prf_error_axis->GetBinCenter(prf_bin+1),
                              sigma);
      _uncertainty_vs_prf_gr->SetPointError(_uncertainty_vs_prf_gr->GetN()-1,
                              0.5*_prf_error_axis->GetBinWidth(prf_bin+1),
                              sigma_e);

    } // loop over prf uncertainty bins
  } // if (_do_separate_fit)

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

  std::cout << "*************** Time consuming **************" << std::endl;
  std::cout << "Reconstruction:\t"        << _sw_partial[0]->CpuTime() * 1.e3 / _EventList.size() << std::endl;
  std::cout << "Analysis:      \t"        << _sw_partial[1]->CpuTime() * 1.e3 / _EventList.size() << std::endl;
  std::cout << "  Col loop:    \t"        << _sw_partial[2]->CpuTime() * 1.e3 / _EventList.size() << std::endl;
  std::cout << "  Fitters:     \t"        << _sw_partial[3]->CpuTime() * 1.e3 / _EventList.size() << std::endl;
  std::cout << "  Filling:     \t"        << _sw_partial[4]->CpuTime() * 1.e3 / _EventList.size() << std::endl;
  return true;
}

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna(argc, argv);
  if (!ana->Initialize())               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}
