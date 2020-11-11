#include "SpatialResolAna.hxx"

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
  _do_separate_pad_fit(false),
  // WARNING
  _correction(false),
  _gaussian_residuals(true),
  _charge_uncertainty(true),
  _gaus_lorentz_PRF(false),
  _iteration(-1){
  //********************************************************************

  // read CLI
  const struct option longopts[] = {
    {"input",           no_argument,    0,    'i'},  // 0
    {"output",          no_argument,    0,    'o'},
    {"batch",           no_argument,    0,    'b'},
    {"verbose",         no_argument,    0,    'v'},
    {"rewrite",         no_argument,    0,    'r'},
    {"correction",      no_argument,    0,    'c'},
    {"full_track_fit",    no_argument,  0,     0},  // 6
    {"separate_pad_fit",  no_argument,  0,     0},  // 7
    {"linear_fit",        no_argument,  0,     0},  // 8
    {"gaus_lorentz",      no_argument,  0,     0}, // 9
    {"help",            no_argument,    0,    'h'},
    {"start",           required_argument,    0,     0},  // 11
    {"end",             required_argument,    0,     0},  // 12
    {0,                 0,              0,     0}
  };

  int index;

  optind = 1;
  for (;;) {
    int c = getopt_long(argc, argv, "i:o:bv:drmst:ca", longopts, &index);
    if (c < 0) break;
    switch (c) {
      case 0 :
        if (index == 6)
          _do_full_track_fit = true;
        if (index == 7)
          _do_separate_pad_fit = true;
        if (index == 8)
          _do_arc_fit = false;
        if (index == 9)
          _gaus_lorentz_PRF = true;
        break;
      case 't' : _iteration        = atoi(optarg);  break;
      case 'c' : _correction       = true;         break;
    }
  }

  if (_iteration == -1) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna().";
    std::cout << " Iteration should be defined as a input param" << std::endl;
    exit(1);
  }

  if (_do_separate_pad_fit && _do_full_track_fit) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna().";
    std::cout << " Incorrect input params definintion" << std::endl;
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

  _uncertainty_vs_prf_gr_prev = NULL;
  _uncertainty_vs_prf_histo   = NULL;

  _PRF_function = InitializePRF("PRF_function");
  _PRF_function_2pad = InitializePRF("PRF_function_2pad");
  _PRF_function_3pad = InitializePRF("PRF_function_3pad");
  _PRF_function_4pad = InitializePRF("PRF_function_4pad");

  if (_do_full_track_fit)
    _PRF_function->FixParameter(0, 1.);

  // load information from previous iteration
  if (_iteration) {
    TString prev_file_name = _file_out_name;
    prev_file_name   = prev_file_name(0, prev_file_name.Index("iter"));
    prev_file_name  += "iter";
    prev_file_name  += TString::Itoa(_iteration - 1, 10);
    prev_file_name  += ".root";

    _Prev_iter_file = new TFile(prev_file_name.Data(), "READ");
    if (!_Prev_iter_file->IsOpen()) {
      std::cerr << "ERROR! SpatialResolAna::Initialize()" << std::endl;
      std::cerr << "File from previous iteration is not found" << std::endl;
      exit(1);
    }
    auto histo_prev = (TH2F*)_Prev_iter_file->Get("PRF_histo");
    histo_prev->SetName("prev_hsto");
    auto gr = new TGraphErrors();
    if (!ProfilePRF(histo_prev, gr)) {
      std::cerr << "ERROR! SpatialResolAna::Initialize()" << std::endl;
      std::cerr << "PRF can not be profiled" << std::endl;
      exit(1);
    }
    gr->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);
    _PRF_function = gr->GetFunction("PRF_function");
    // _PRF_function   = (TF1*)_Prev_iter_file->Get("PRF_function");
    auto uncertainty_graph = (TH1F*)_Prev_iter_file->Get("resol_total");

    if (!_PRF_function || !uncertainty_graph) {
      std::cerr << "ERROR. SpatialResolAna::Initialize().";
      std::cout << "PRF function or resolution is not specified" << std::endl;
      std::cerr << "Search in " << prev_file_name << std::endl;
      exit(1);
    }

    // if (_PRF_function) {
    //   std::cout << "      PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula() << "  with ";
    //   for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
    //     std::cout << "  " << _PRF_function->GetParameter(i) << ",";
    //   std::cout << "  Chi2/NDF " << _PRF_function->GetChisquare()
    //             << "/" << _PRF_function->GetNDF() << std::endl;
    //   std::cout << std::endl;
    // }

    Double_t mean, sigma;
    sigma = 0.5 * GetFWHM(uncertainty_graph, mean);
    uncertainty_graph->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);
    _uncertainty = uncertainty_graph->GetFunction("gaus")->GetParameter(2);

    TGraphErrors* temp = NULL;
    if (_do_separate_pad_fit)
      temp = (TGraphErrors*)_Prev_iter_file->Get("uncertainty_vs_prf_gr");
    // transform graph into histo
    if (temp && temp->GetN()) {
      _uncertainty_vs_prf_gr_prev = temp;
      _uncertainty_vs_prf_histo = new TH1F("uncertainty_histo", "",
          prf_error_bins-1, prf_error_bins_arr);
      for (auto bin_id = 1; bin_id <= _prf_error_axis->GetNbins(); ++bin_id) {
        _uncertainty_vs_prf_histo->SetBinContent(
          bin_id,
          _uncertainty_vs_prf_gr_prev->GetY()[bin_id-1]);
      }
    }

    // read event list passed through reconstruction+selection
    // at the moment skip for the 1 iteration (after 0) files with TEvent
    // reason: in at the step 0 TEvent file is generated and the list of events
    // is different from the initial (raw) data file
    // need to look through all events again
    if (!_work_with_event_file) {
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
  } // if iteration

  _qulity_ratio   = new TH1F("quality_ratio",
    "Quality ratio arc/linear", 100, 0., 2.);
  _chi2_ratio     = new TH1F("chi_ratio",
    "Chi2 ratio arc/linear", 100, 0., 2.);

  _mom_reco   = new TH1F("mom_reco", "", 2000, -6., 6.);
  _pos_reco   = new TH1F("pos_reco", "", 8000, -0.2, 0.2);
  _ang_reco   = new TH1F("ang_reco", "", 3000, -0.3, 0.3);

  // Initialise histoes and graphs
  _PRF_histo = new TH2F("PRF_histo","", prf_bin, prf_min, prf_max, 150,0.,1.5);

  // Initialize graph for PRF profiling
  _PRF_graph = new TGraphErrors();
  _PRF_graph->SetName("PRF_graph");
  _PRF_graph_2pad = new TGraphErrors();
  _PRF_graph_2pad->SetName("PRF_graph_2pad");
  _PRF_graph_3pad = new TGraphErrors();
  _PRF_graph_3pad->SetName("PRF_graph_3pad");
  _PRF_graph_4pad = new TGraphErrors();
  _PRF_graph_4pad->SetName("PRF_graph_4pad");

  // PRF for different multiplicities
  auto dir_prf_mult = _file_out->mkdir("prf_mult");
  _output_vector.push_back(dir_prf_mult);

  _file_out->cd();
  _tree = new TTree("outtree", "");
  _tree->Branch("ev",           &_ev);
  _tree->Branch("angle_yz",     &_angle_yz);
  _tree->Branch("angle_xy",     &_angle_xy);
  _tree->Branch("multiplicity",
                &_multiplicity,
                TString::Format("multiplicity[%i]/I", geom::nPadx)
                );
  _tree->Branch("charge",
                &_charge,
                TString::Format("charge[%i]/I", geom::nPadx)
                );
  _tree->Branch("residual",
                &_residual,
                TString::Format("residual[%i]/F", geom::nPadx)
                );
  _tree->Branch("dx",
                &_dx,
                TString::Format("dx[%i][10]/F", geom::nPadx)
                );
  _tree->Branch("qfrac",
                &_qfrac,
                TString::Format("qfrac[%i][10]/F", geom::nPadx)
                );
  _tree->Branch("clust_pos",
                &_clust_pos,
                TString::Format("clust_pos[%i]/I", geom::nPadx)
                );
  _tree->Branch("track_pos",
                &_track_pos,
                TString::Format("track_pos[%i]/I", geom::nPadx)
                );

  _output_vector.push_back(_tree);

  _PRF_histo_2pad = new TH2F("PRF_histo_2pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _PRF_histo_3pad = new TH2F("PRF_histo_3pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _PRF_histo_4pad = new TH2F("PRF_histo_4pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);

  dir_prf_mult->Append(_PRF_histo_2pad);
  dir_prf_mult->Append(_PRF_histo_3pad);
  dir_prf_mult->Append(_PRF_histo_4pad);

  dir_prf_mult->Append(_PRF_function_2pad);
  dir_prf_mult->Append(_PRF_function_3pad);
  dir_prf_mult->Append(_PRF_function_4pad);

  dir_prf_mult->Append(_PRF_graph_2pad);
  dir_prf_mult->Append(_PRF_graph_3pad);
  dir_prf_mult->Append(_PRF_graph_4pad);

  for (auto i = 0; i < geom::GetMaxColumn(_invert); ++i)
    _PRF_histo_col[i] = new TH2F(Form("PRF_histo_col_%i", i),
      "", prf_bin, prf_min, prf_max, 150,0.,1.5);

  for (auto i = 0; i < 4; ++i) {
    _PRF_histo_xscan[i] = new TH2F(Form("PRF_histo_pad_%i", i),
      "", prf_bin, prf_min, prf_max, 150,0.,1.5);
    _PRF_graph_xscan[i] = new TGraphErrors();
    _PRF_graph_xscan[i]->SetName(Form("PRF_graph_pad_%i", i));
  }

  _resol_total = new TH1F("resol_total",
     "", resol_bin, resol_min, resol_max);

  _output_vector.push_back(_resol_total);

  for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
    _resol_col_hist[j]  = new TH1F(Form("resol_histo_%i", j),
     "", resol_bin, resol_min, resol_max);
    _resol_col_hist_except[j]  = new TH1F(Form("resol_histo1__%i", j),
     "", resol_bin, resol_min, resol_max);

    _resol_col_hist_2pad[j]  = new TH1F(Form("resol_histo_2pad_%i", j),
      "", resol_bin, resol_min, resol_max);
    _resol_col_hist_2pad_except[j]  = new TH1F(Form("resol_histo1_2pad_%i", j),
     "", resol_bin, resol_min, resol_max);

    _resol_col_hist_3pad[j]  = new TH1F(Form("resol_histo_3pad_%i", j),
      "", resol_bin, resol_min, resol_max);
    _resol_col_hist_3pad_except[j]  = new TH1F(Form("resol_histo1_3pad_%i", j),
     "", resol_bin, resol_min, resol_max);
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

  _uncertainty_vs_prf_gr = new TGraphErrors();
  _uncertainty_vs_prf_gr->SetName("uncertainty_vs_prf_gr");
  _output_vector.push_back(_uncertainty_vs_prf_gr);

  auto dir_resol = _file_out->mkdir("resol_column");
  _output_vector.push_back(dir_resol);
  for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
    dir_resol->Append(_resol_col_hist[j]);
    dir_resol->Append(_resol_col_hist_except[j]);
  }

  if (_do_separate_pad_fit) {
    Double_t arr[4] = {0., .15, .7, 1.01};
    _prf_scale_axis = new TAxis(3, arr);
    for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
      for (auto id=0; id < 3; ++id) {
        _Fit_quality_plots[id][j] = new TH1F(Form("pad_fit_q_%i_%i", j, id),
         "", 300, 0., 30.);
        _output_vector.push_back(_Fit_quality_plots[id][j]);
      }
    }
  }

  auto dir_x_scan = _file_out->mkdir("x_scan");
  _output_vector.push_back(dir_x_scan);

  _x_scan_axis = new TAxis(x_scan_bin, x_scan_min, x_scan_max);
  for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
    for (auto i = 0; i < x_scan_bin; ++i) {
      _resol_col_x_scan[j][i] = new TH1F(Form("resol_histo_Xscan_%i_%i", j, i),
       "", resol_bin, resol_min, resol_max);
      dir_x_scan->Append(_resol_col_x_scan[j][i]);
      _mult_x_scan[j][i] = new TH1F(Form("mult_histo_Xscan_%i_%i", j, i),
       "multiplicity", 10, 0., 10.);
      dir_x_scan->Append(_mult_x_scan[j][i]);
    }
  }

  for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
    for (auto i = 0; i < x_scan_bin; ++i) {
      _resol_col_x_scan_lim_mult[j][i] =
        new TH1F(Form("resol_histo_Xscan_%i_%i", j, i), "",
        resol_bin, resol_min, resol_max);
      dir_x_scan->Append(_resol_col_x_scan_lim_mult[j][i]);
    }
  }

  for (auto i = 0; i < prf_error_bins-1; ++i) {
    _uncertainty_prf_bins[i] = new TH1F(Form("error_prf_bin_%i", i),
      "", resol_bin, resol_min, resol_max);
    // _output_vector.push_back(_uncertainty_prf_bins[i]);
  }

  auto dir_prf = _file_out->mkdir("prf_column");
  _output_vector.push_back(dir_prf);
  for (auto j = 0; j < geom::GetMaxColumn(_invert); ++j) {
    dir_prf->Append(_PRF_histo_col[j]);
  }

  _passed_events.clear();

  std::cout << "done" << std::endl;
  if (_verbose > 1) {
    std::cout << "\t PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula();
    std::cout << "  with ";
    for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
      std::cout << "  " << _PRF_function->GetParameter(i) << ",";
    std::cout << std::endl;
  }

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  // Initialise track fitter
  _fitter = new TrackFitter(TrackFitter::CERN_like, _PRF_function,
      fit_bound_right, _uncertainty, _iteration, _verbose, _invert,
      _charge_uncertainty, _do_arc_fit);

  if (_do_full_track_fit)
    _fitter->SetType(TrackFitter::ILC_like);
  else if (_do_separate_pad_fit)
    _fitter->SetType(TrackFitter::Separate_pads);

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
bool SpatialResolAna::ProcessEvent(const TEvent* event) {
  //********************************************************************

  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    if (_verbose > 1)
      std::cout << "Track id = " << trackId << std::endl;

    if (!sel::CrossingTrackSelection(track, _invert, _verbose))
      continue;

    std::vector<double> fit_v = sel::GetFitParams(track, _invert);
    std::vector<double> fit_xz = sel::GetFitParamsXZ(track, _invert);

    _angle_xy = fit_v[2];
    _angle_yz = fit_xz[2] * sel::v_drift_est;

    _ev = event->GetID();

    _store_event = true;

    int cluster[geom::nPadx];
    int cluster_N[geom::nPadx];
    double track_pos[geom::nPadx];
    double cluster_mean[geom::nPadx];
    float charge_max[geom::nPadx];
    double a_peak_fit[geom::nPadx];

    pads_t pos_in_pad;
    pos_in_pad.clear();
    pos_in_pad.resize(geom::GetMaxColumn(_invert));
    int Ndots = 0;

    // At the moment ommit first and last column
    // first loop over columns
    _sw_partial[2]->Start(false);

    // reset tree values
    for (auto colId = 0; colId < geom::nPadx; ++colId) {
      _multiplicity[colId]  = -999;
      _charge[colId]        = -999;
      _residual[colId]      = -999;
      for (auto padId = 0; padId < 10; ++padId) {
        _dx[colId][padId]   = -999;
        _qfrac[colId][padId] = -999;
      }

      cluster[colId]       = 0;
      cluster_N[colId]     = 0;
      charge_max[colId]    = 0;
      track_pos[colId]     = -999.;
      cluster_mean[colId]  = 0.;
      a_peak_fit[colId]    = 0.;
    }

    auto robust_cols = GetRobustCols(track->GetCols(_invert));
    for (auto col:robust_cols) {
      if (!col[0])
        continue;
      auto it_x = col[0]->GetCol(_invert);

      TH1F* cluster_h;
      if (!_invert)
        cluster_h = new TH1F("cluster", "", geom::nPady,
          -1.*geom::nPady*geom::dy/2., 1.*geom::nPady*geom::dy/2.);
      else
        cluster_h = new TH1F("cluster", "", geom::nPadx,
          -1*geom::nPadx*geom::dx/2., 1.*geom::nPadx*geom::dx/2.);

      // loop over rows
      auto robust_pads = GetRobustPadsInColumn(col);
      _multiplicity[it_x] = robust_pads.size();
      _charge[it_x] = std::accumulate(robust_pads.begin(), robust_pads.end(), 0,
                                      [](const int& x, const THit* hit)
                                      {return x + hit->GetQ();});

      for (auto pad:robust_pads) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow(_invert);
        auto q = pad->GetQ();

        cluster[it_x] += q;
        ++cluster_N[it_x];
        cluster_h->Fill(geom::GetYpos(it_y, _invert), q);
        if (charge_max[it_x] < q)
          charge_max[it_x] = q;
      } // end of loop over rows

      cluster_mean[it_x] = cluster_h->GetMean();

      delete cluster_h;
      cluster_h = NULL;

      if (!cluster[it_x] || !charge_max[it_x])
        continue;

      ++Ndots;

      track_pos[it_x] = _fitter->FitCluster(robust_pads,
          cluster[it_x], cluster_mean[it_x], pos_in_pad,
          _uncertainty_vs_prf_histo);

      if (_verbose > 5) {
        std::cout << "Cluster pos " << track_pos[it_x] << "\t";
        std::cout << cluster_mean[it_x] << std::endl;
      }
    } // loop over columns

    if (_verbose > 3)
      std::cout << "Loop over columns done" << std::endl;

    _Cols_used->Fill(Ndots);

    _sw_partial[2]->Stop();
    _sw_partial[3]->Start(false);

    TF1* fit = NULL;
    fit = _fitter->FitTrack(track_pos, cluster_N, track,
                            track_pos[1], pos_in_pad);

    if (!fit)
      continue;

    // TODO review this mess
    if (fit && _do_full_track_fit)
      fit = (TF1*)fit->Clone();

    double quality = fit->GetChisquare() / fit->GetNDF();
    _Chi2_track->Fill(quality);

    if (_verbose > 3)
      std::cout << "Track fit done" << std::endl;

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
    for (int i = 1; i < geom::GetMaxColumn(_invert)-1; ++i) {
      if (!_correction)
        fit1[i] = fit;
      else {
        fit1[i] = _fitter->FitTrack(track_pos, cluster_N, track, track_pos[1],
          pos_in_pad, i);
        // TODO review this mess
        if (_do_full_track_fit)
          fit1[i] = (TF1*)fit1[i]->Clone();
      }
    }

    _sw_partial[3]->Stop();
    _sw_partial[4]->Start(false);

    // second loop over columns
    for (auto col:robust_cols) {
      if (!col[0])
        continue;
      auto it_x = col[0]->GetCol(_invert);

      if (track_pos[it_x]  == -999.)
        continue;

      double x    = geom::GetXpos(it_x, _invert);
      double track_fit_y    = fit->Eval(x);
      double track_fit_y1   = fit1[it_x]->Eval(x);

      if (_verbose > 4) {
        std::cout << "Residuals cluster:track\t" << track_pos[it_x] << "\t";
        std::cout << track_fit_y << "\t" << (track_pos[it_x] - track_fit_y)*1e6;
        std::cout << std::endl;
      }

      // fill SR
      _clust_pos[it_x] = track_pos[it_x];
      _track_pos[it_x] = track_fit_y;
      _residual[it_x] = track_pos[it_x] - track_fit_y;
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

      // if (cluster_N[it_x] == 2) {
      //   _resol_col_hist[it_x]->Fill(track_pos[it_x] - track_fit_y);
      //   _resol_col_hist_2pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      // } else if (cluster_N[it_x] == 3) {
      //   _resol_col_hist_3pad[it_x]->Fill(track_pos[it_x] - track_fit_y);
      //   _resol_col_hist_3pad_except[it_x]->Fill(track_pos[it_x] - track_fit_y1);
      // }

      // in case of full track fit calc normalisation coefficient
      // if (_do_full_track_fit) {
      //   float a_nom = 0.;
      //   float a_den = 0.;
      //   for (auto pad:col) {
      //     auto q      = pad->GetQ();
      //     auto it_y   = pad->GetRow(_invert);
      //     if (!q)
      //       continue;
      //     double center_pad_y = geom::GetYpos(it_y, _invert);

      //     if (abs(center_pad_y - track_fit_y) > fit_bound_right)
      //       continue;

      //     a_nom += _PRF_function->Eval(center_pad_y - track_fit_y);
      //     a_den += TMath::Power(_PRF_function->Eval(center_pad_y - track_fit_y), 2) /   q;
      //   }

      //   a_peak_fit[it_x] = a_nom / a_den;
      // } else
      a_peak_fit[it_x] = 1.*cluster[it_x];

      // fill pad accuracy
      // if (_do_separate_pad_fit) {
      //   for (auto it = pos_in_pad[it_x].begin();
      //             it < pos_in_pad[it_x].end();
      //             ++it) {
      //     int bin_prf = -1 + _prf_scale_axis->FindBin((*it).first);
      //   // TODO make the definition through constants
      //     if (bin_prf < 0 || bin_prf > 2) {
      //       std::cout << "Error bin 1  " << bin_prf << "\t" << (*it).first << std::endl;
      //       exit(1);
      //     }
      //     _Fit_quality_plots[bin_prf][it_x]->Fill(abs(track_fit_y - (*it).second.first) / (*it).second.second);

      //     int bin2 = -1 + _prf_error_axis->FindBin((*it).first);
      //     if (bin2 < 0 || bin2 >= prf_error_bins-1) {
      //       std::cout << "Error bin 2  " << bin2 << "\t" << (*it).first << std::endl;
      //       exit(1);
      //     }
      //     _uncertainty_prf_bins[bin2]->Fill(track_fit_y - (*it).second.first);
      //   }
      // }

      if (cluster_N[it_x] == 1)
        continue;
      // Fill PRF
      auto robust_pads = GetRobustPadsInColumn(col);
      int padId = 0;
      for (auto pad:robust_pads) {
        if (!pad)
          continue;

        auto it_y = pad->GetRow(_invert);
        auto q = pad->GetQ();

        if (!cluster[it_x] || !q)
          continue;

        double center_pad_y = geom::GetYpos(it_y, _invert);

        if (it_x == 4) {
          auto bin_prf = _x_pads->FindBin(center_pad_y);
          if (bin_prf > 0 && bin_prf < 5) {
            _PRF_histo_xscan[bin_prf-1]->Fill(center_pad_y - track_fit_y,
                                              q / a_peak_fit[it_x]);
          }
        }

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

        if (padId > 9)
          continue;

        _dx[it_x][padId]    = center_pad_y - track_fit_y;
        _qfrac[it_x][padId] = q / a_peak_fit[it_x];

        // robust_pads are assumed sorted!!
        ++padId;
      }
    } // loop over colums
    _sw_partial[4]->Stop();
    for (int i = 0; i < geom::GetMaxColumn(_invert); ++i)
      if (fit1[i] && _correction) {
        delete fit1[i];
        fit1[i] = NULL;
      }

    delete fit;

    if(_test_mode) this->DrawSelectionCan(event,trackId);
  } // loop over tracks

  _tree->Fill();

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

  std::cout << "PRF profiling............................";

  ProfilePRF(_PRF_histo, _PRF_graph);
  ProfilePRF(_PRF_histo_2pad, _PRF_graph_2pad);
  ProfilePRF(_PRF_histo_3pad, _PRF_graph_3pad);
  ProfilePRF(_PRF_histo_4pad, _PRF_graph_4pad);

  std::cout << "done" << std::endl;
  std::cout << "PRF fit..................................";

  _PRF_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);
  _PRF_graph_2pad->Fit("PRF_function_2pad", "Q", "", -0.014, 0.014);
  _PRF_graph_3pad->Fit("PRF_function_3pad", "Q", "", fit_bound_left, fit_bound_right);
  _PRF_graph_4pad->Fit("PRF_function_4pad", "Q", "", fit_bound_left, fit_bound_right);

  std::cout << "done" << std::endl;
  std::cout << "Process x histoes........................";

  for (auto i = 0; i < 4; ++i)
    ProfilePRF(_PRF_histo_xscan[i], _PRF_graph_xscan[i]);

  std::cout << "done" << std::endl;
  std::cout << "Process histoes..........................";

  TH1F* resol = new TH1F("resol", "", 1000, 0., 0.001);

  for (auto i = 0; i <= geom::GetMaxColumn(_invert) - 1; ++i)
      _resol_total->Add(_resol_col_hist[i]);

  for (auto i = 1; i < geom::GetMaxColumn(_invert) - 1; ++i) {

    TH1F* res     = _resol_col_hist[i];
    TH1F* res_e   =  _resol_col_hist_except[i];

    Double_t mean, sigma, sigma_ex;
    Double_t mean_e = 0., sigma_e = 0., sigma_ex_e = 0.;

    sigma = 0.5 * GetFWHM(res, mean);

    if (!res->Integral()) {
      std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Empty residuals" << std::endl;
      exit(1);
    }

    if (_gaussian_residuals) {
      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);
      _resol_col_hist_2pad[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad[i]->Fit("gaus", "Q");



      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        continue;
      }

      res_e->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;
      _resol_col_hist_2pad_except[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad_except[i]->Fit("gaus", "Q");

      TF1* func_ex = res_e->GetFunction("gaus");
      if (!func_ex) {
        std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Exeptional residual fit fail" << std::endl;
        continue;
      }

      mean     = func->GetParameter(1);
      sigma    = func->GetParameter(2);
      sigma_ex = func_ex->GetParameter(2);

      mean_e      = func->GetParError(1);
      sigma_e     = func->GetParError(2);
      sigma_ex_e  = func_ex->GetParError(2);
    } else {
      // use FWHM
      sigma_ex = 0.5 * GetFWHM(res_e, mean);
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

  if (_do_separate_pad_fit && _iteration) {
    for (auto prf_bin = 0; prf_bin < prf_error_bins-1; ++prf_bin) {
      TH1F* res = _uncertainty_prf_bins[prf_bin];

      Double_t mean;
      Double_t sigma = 0.5*GetFWHM(res, mean);

      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;

      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cerr << "ERROR. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        exit(1);
      }

      Double_t sigma_e;

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

  if (_PRF_graph->GetFunction("PRF_function")) {
    _PRF_function = _PRF_graph->GetFunction("PRF_function");
    std::cout << "      PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula() << "  with ";
    for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
      std::cout << "  " << _PRF_function->GetParameter(i) << ",";
    std::cout << "  Chi2/NDF " << _PRF_graph->GetFunction("PRF_function")->GetChisquare()
              << "/" << _PRF_graph->GetFunction("PRF_function")->GetNDF() << std::endl;
    std::cout << std::endl;

    _resol_total->Fit("gaus", "Q");
    auto func = _resol_total->GetFunction("gaus");
    TString output = "NaN";
    if (func)
      output = TString().Itoa(1.e6*func->GetParameter(2), 10) + " um";

    std::cout << "Spatial resolution\t" << output << std::endl;
  }

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

bool SpatialResolAna::ProfilePRF(const TH2F* PRF_h, TGraphErrors* gr) {
  if (!PRF_h)
    return false;

  // Float_t threshold = 0.05 * PRF_h->GetMaximum();
  Float_t threshold = 0.;

  for (auto i = 1; i < PRF_h->GetXaxis()->GetNbins(); ++i) {

    TH1D* temp_h = PRF_h->ProjectionY(Form("projections_bin_%i", i), i, i);
    if (!temp_h)
      return false;

    if (temp_h->GetMaximum() < threshold)
      continue;

    double x = PRF_h->GetXaxis()->GetBinCenter(i);
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

    gr->SetPoint(gr->GetN(), x, y);
    gr->SetPointError(gr->GetN()-1, 0, e/2.);
  } // end of PRF histo profiling

  return true;
}

Double_t SpatialResolAna::GetFWHM(const TH1F* h, Double_t& mean) {
  mean = h->GetMean();
  float max   = h->GetMaximum();
  float start = h->GetBinLowEdge(h->FindFirstBinAbove(max/2));
  float end   = h->GetBinLowEdge(h->FindLastBinAbove(max/2)) +
  h->GetBinWidth(h->FindLastBinAbove(max/2));

  return end - start;
}

TF1* SpatialResolAna::InitializePRF(const TString name) {
  TF1* func;
  if (_gaus_lorentz_PRF) {
    func = new TF1(name,
      "[0] * exp(-4*TMath::Log(2)*(1-[1])*TMath::Power(x/[2], 2.)) / (1+4 * [1] * TMath::Power(x/[2], 2.) )",
      prf_min, prf_max);
    func->SetParName(0, "Const");
    func->SetParName(1, "r");
    func->SetParName(2, "w");

    auto c = 1.;
    auto r = 0.5;
    auto s = 0.005;
    func->SetParameters(c, r, s);

    return func;
  }

  func  = new TF1(name,
    "[0]*(1+[1]*x*x + [2] * x*x*x*x) / (1+[3]*x*x+[4]*x*x*x*x)",
    prf_min, prf_max);
  func->SetParName(0, "Const");
  func->SetParName(1, "a2");
  func->SetParName(2, "a4");
  func->SetParName(3, "b2");
  func->SetParName(4, "b4");

  double co = 1.;
  double a2 = 2.35167e3;
  double a4 = 6.78962e7;
  double b2 = 3.36748e3;
  double b4 = 6.45311e8;
  func->SetParameters(co, a2, a4, b2, b4);

  return func;
}

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna(argc, argv);
  if (!ana->Initialize())               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}
