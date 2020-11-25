#include "SpatialResolAna.hxx"

// TODO enum?
//! verbosity level
//! 1 (default) print progress, memory usage, time consumption
//! 2           print event and track number
//! 3           print analysis steps
//! 4           print fit details
//! 5           print residuals
//! 6           print PRF details

//******************************************************************************
SpatialResolAna::SpatialResolAna(int argc, char** argv):
  AnalysisBase(argc, argv),
  _Prev_iter_name(TString("")),
  _do_linear_fit(false),
  _do_para_fit(false),
  _do_full_track_fit(false),
  _do_separate_pad_fit(false),
  // WARNING
  _correction(false),
  _gaussian_residuals(true),
  _charge_uncertainty(true),
  _gaus_lorentz_PRF(false),
  _diagonal(false),
  _iteration(-1){
//******************************************************************************

  // read CLI
  const struct option longopts[] = {
    {"input",           no_argument,    0,    'i'},         // 0
    {"output",          no_argument,    0,    'o'},         // 1
    {"batch",           no_argument,    0,    'b'},         // 2
    {"verbose",         no_argument,    0,    'v'},         // 3
    {"rewrite",         no_argument,    0,    'r'},         // 4
    {"correction",      no_argument,    0,    'c'},         // 5
    {"start",           required_argument,      0,     0},  // 6
    {"end",             required_argument,      0,     0},  // 7
    // fitter method
    {"full_track_fit",  no_argument,    0,      0},         // 8
    {"separate_pad_fit",no_argument,    0,      0},         // 9
    // track shape
    {"linear_fit",      no_argument,    0,      0},         // 10
    {"para_fit",        no_argument,    0,      0},         // 11
    // PRF  shape
    {"gaus_lorentz",    no_argument,    0,      0},         // 12

    {"diagonal",        no_argument,    0,      0},         // 13

    // {"previous",        required_argument, 0, 'p'},         // 14
    {"help",            no_argument,    0,    'h'},         // 15

    {0,                 0,              0,      0}
  };

  int index;
  optind = 1;
  for (;;) {
    int c = getopt_long(argc, argv, "i:o:bv:drmst:cap:", longopts, &index);
    if (c < 0) break;
    switch (c) {
      case 0 :
        if (index == 8)
          _do_full_track_fit = true;
        if (index == 9)
          _do_separate_pad_fit = true;
        if (index == 10)
          _do_linear_fit = true;
        if (index == 11)
          _do_para_fit = true;
        if (index == 12)
          _gaus_lorentz_PRF = true;
        if(index == 13)
          _diagonal = true;
        break;
      case 't' : _iteration        = atoi(optarg);  break;
      case 'c' : _correction       = true;          break;
      case 'p' :  _Prev_iter_name  = optarg;        break;
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

//******************************************************************************
bool SpatialResolAna::Initialize() {
//******************************************************************************
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
  // _PRF_function_2pad = InitializePRF("PRF_function_2pad");
  // _PRF_function_3pad = InitializePRF("PRF_function_3pad");
  // _PRF_function_4pad = InitializePRF("PRF_function_4pad");

  if (_do_full_track_fit)
    _PRF_function->FixParameter(0, 1.);

  // load information from previous iteration
  if (_iteration) {
    if (_Prev_iter_name.Length() == 0) {
      _Prev_iter_name = _file_out_name;
      _Prev_iter_name   = _Prev_iter_name(0, _Prev_iter_name.Index("iter"));
      _Prev_iter_name  += "iter";
      _Prev_iter_name  += TString::Itoa(_iteration - 1, 10);
      _Prev_iter_name  += ".root";
    }

    _Prev_iter_file = new TFile(_Prev_iter_name.Data(), "READ");
    if (!_Prev_iter_file->IsOpen()) {
      std::cerr << "ERROR! SpatialResolAna::Initialize()" << std::endl;
      std::cerr << "File from previous iteration is not found" << std::endl;
      std::cerr << "File name: " << _Prev_iter_name << std::endl;
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
    auto uncertainty_graph = (TH1F*)_Prev_iter_file->Get("resol_total");

    if (!_PRF_function || !uncertainty_graph) {
      std::cerr << "ERROR. SpatialResolAna::Initialize().";
      std::cout << "PRF function or resolution is not specified" << std::endl;
      std::cerr << "Search in " << _Prev_iter_name << std::endl;
      exit(1);
    }

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
                TString::Format("multiplicity[%i]/I", Nclusters)
                );
  _tree->Branch("charge",
                &_charge,
                TString::Format("charge[%i]/I", Nclusters)
                );
  _tree->Branch("residual",
                &_residual,
                TString::Format("residual[%i]/F", Nclusters)
                );
  _tree->Branch("dx",
                &_dx,
                TString::Format("dx[%i][10]/F", Nclusters)
                );
  _tree->Branch("qfrac",
                &_qfrac,
                TString::Format("qfrac[%i][10]/F", Nclusters)
                );
  _tree->Branch("time",
                &_time,
                TString::Format("time[%i][10]/I", Nclusters)
                );
  _tree->Branch("clust_pos",
                &_clust_pos,
                TString::Format("clust_pos[%i]/F", Nclusters)
                );
  _tree->Branch("track_pos",
                &_track_pos,
                TString::Format("track_pos[%i]/F", Nclusters)
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

  // dir_prf_mult->Append(_PRF_function_2pad);
  // dir_prf_mult->Append(_PRF_function_3pad);
  // dir_prf_mult->Append(_PRF_function_4pad);

  dir_prf_mult->Append(_PRF_graph_2pad);
  dir_prf_mult->Append(_PRF_graph_3pad);
  dir_prf_mult->Append(_PRF_graph_4pad);

  for (auto i = 0; i < Nclusters; ++i)
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

  for (auto j = 0; j < Nclusters; ++j) {
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
  for (auto j = 0; j < Nclusters; ++j) {
    dir_resol->Append(_resol_col_hist[j]);
    dir_resol->Append(_resol_col_hist_except[j]);
  }

  if (_do_separate_pad_fit) {
    Double_t arr[4] = {0., .15, .7, 1.01};
    _prf_scale_axis = new TAxis(3, arr);
    for (auto j = 0; j < Nclusters; ++j) {
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
  for (auto j = 0; j < Nclusters; ++j) {
    for (auto i = 0; i < x_scan_bin; ++i) {
      _resol_col_x_scan[j][i] = new TH1F(Form("resol_histo_Xscan_%i_%i", j, i),
       "", resol_bin, resol_min, resol_max);
      dir_x_scan->Append(_resol_col_x_scan[j][i]);
      _mult_x_scan[j][i] = new TH1F(Form("mult_histo_Xscan_%i_%i", j, i),
       "multiplicity", 10, 0., 10.);
      dir_x_scan->Append(_mult_x_scan[j][i]);
    }
  }

  for (auto j = 0; j < Nclusters; ++j) {
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
  for (auto j = 0; j < Nclusters; ++j) {
    dir_prf->Append(_PRF_histo_col[j]);
  }

  _passed_events.clear();

  std::cout << "done" << std::endl;
  if (_verbose >= v_event_number) {
    std::cout << "\t PRF(x) = " << _PRF_function->GetFormula()->GetExpFormula();
    std::cout << "  with ";
    for (auto i = 0; i < _PRF_function->GetNpar(); ++i)
      std::cout << "  " << _PRF_function->GetParameter(i) << ",";
    std::cout << std::endl;
  }

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize(_verbose);

  // Initialise track fitter
  TrackFitterBase::TrackShape shape = TrackFitterBase::arc;
  if (_do_linear_fit) {
    shape = TrackFitterBase::linear;
  } else if (_do_para_fit) {
    shape = TrackFitterBase::parabola;
  }

  _fitter = new TrackFitCern(shape,
                             _invert,
                             _diagonal,
                             _verbose,
                             _iteration,
                             _PRF_function,
                             fit_bound_right,
                             _charge_uncertainty
                             );

  // Initialize timers
  _sw_partial[2] = new TStopwatch();
  _sw_partial[2]->Reset();
  _sw_partial[3] = new TStopwatch();
  _sw_partial[3]->Reset();
  _sw_partial[4] = new TStopwatch();
  _sw_partial[4]->Reset();

  return true;
}

//******************************************************************************
bool SpatialResolAna::ProcessEvent(const TEvent* event) {
//******************************************************************************

  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    if (_verbose >= v_event_number)
      std::cout << "Track id = " << trackId << std::endl;

    if (!sel::CrossingTrackSelection(track, _invert, _verbose))
      continue;

    std::vector<double> fit_v = sel::GetFitParams(track, _invert);
    std::vector<double> fit_xz = sel::GetFitParamsXZ(track, _invert);

    _angle_xy = fit_v[2];
    _angle_yz = fit_xz[2] * sel::v_drift_est;

    _ev = event->GetID();

    _store_event = true;

    int cluster[Nclusters];
    int cluster_N[Nclusters];
    double track_pos[Nclusters];
    double cluster_mean[Nclusters];
    float charge_max[Nclusters];
    double a_peak_fit[Nclusters];

    pads_t pos_in_pad;
    pos_in_pad.clear();
    pos_in_pad.resize(Nclusters);
    int Ndots = 0;

    // At the moment ommit first and last column
    // first loop over columns
    _sw_partial[2]->Start(false);

    // reset tree values
    for (auto colId = 0; colId < Nclusters; ++colId) {
      _multiplicity[colId]  = -999;
      _charge[colId]        = -999;
      _residual[colId]      = -999;
      _clust_pos[colId]     = -999;
      _track_pos[colId]     = -999;
      _x[colId]             = -999;
      _x_av[colId]          = -999;
      _cluster_av[colId]    = -999;
      for (auto padId = 0; padId < 10; ++padId) {
        _dx[colId][padId]   = -999;
        _time[colId][padId]   = -999;
        _qfrac[colId][padId] = -999;
      }

      cluster[colId]       = 0;
      cluster_N[colId]     = 0;
      charge_max[colId]    = 0;
      track_pos[colId]     = -999.;
      cluster_mean[colId]  = 0.;
      a_peak_fit[colId]    = 0.;
    }
// *****************************************************************************
// ************** See steps documentation in the README.md file ****************
// *******************  STEP 1 *************************************************

    std::vector<TCluster*> clusters;
    if (_diagonal)
      clusters = DiagonolizeTrack(track);
    else
      clusters = ColonizeTrack(track);

    if (_verbose >= v_analysis_steps)
      std::cout << "Clusterization done" << std::endl;

    if (_diagonal) {
      if (clusters.size() < 10)
        continue;
      // clean first and last column
      sort(clusters.begin(), clusters.end(),
           [&](TCluster* cl1, TCluster* cl2){
              return cl1->GetX() < cl2->GetX();
            });
      clusters.erase(clusters.begin());
      clusters.erase(clusters.end()-1);
    }
    // truncation
    auto robust_clusters = GetRobustCols(clusters);

    if (_verbose >= v_analysis_steps)
      std::cout << "clearing done, columns\t" << robust_clusters.size() << std::endl;

    if (_diagonal && robust_clusters.size() < 40)
      continue;
// *******************  STEP 2 *************************************************

    if (_verbose >= v_analysis_steps)
      std::cout << "start cluster fit" << std::endl;

    for (uint clusterId = 0; clusterId < robust_clusters.size(); ++clusterId) {
      auto cluster = robust_clusters[clusterId];
      if (!cluster->GetHits()[0])
        continue;

      // loop over rows
      auto robust_pads = GetRobustPadsInColumn(cluster->GetHits());
      _multiplicity[clusterId] = robust_pads.size();

      _clust_pos[clusterId] = 0.;
      _charge[clusterId] = 0;
      for (auto pad:robust_pads) {
        if (!pad)
          continue;

        _clust_pos[clusterId] += pad->GetQ() * geom::GetYposPad(pad,
                                                          _invert,
                                                          _diagonal ? units::a45 : 0);
        _charge[clusterId] += pad->GetQ();
      } // loop over pads

      _clust_pos[clusterId] /= _charge[clusterId];
      _x[clusterId] = geom::GetXposPad(robust_pads[0],
                                       _invert,
                                       _diagonal ? units::a45 : 0);

      // WARNING tmp
      double CoC =  _clust_pos[clusterId];

      ++Ndots;

      if (_multiplicity[clusterId] > 1 && _iteration > 0){
        track_pos[clusterId] = _fitter->FitCluster(robust_pads,
                                                   _charge[clusterId],
                                                   _clust_pos[clusterId]
                                                   );
        _clust_pos[clusterId] = track_pos[clusterId];
      } else {
        track_pos[clusterId] = _clust_pos[clusterId];
      }

      cluster->SetY(track_pos[clusterId]);
      cluster->SetCharge(_charge[clusterId]);

      if (_verbose >= v_fit_details) {
        for (auto pad:cluster->GetHits()) {
          std::cout << pad->GetRow(_invert) <<  " : " << pad->GetCol(_invert) << "\t";
        }
        std::cout << std::endl;
      }

      if (_verbose >= v_fit_details)
        std::cout << "X:CoC:Fit\t" << cluster->GetX() << "\t" << CoC << "\t" << cluster->GetY() << std::endl;


      // TODO review this mess
      if (cluster->GetHits().size() == 1) {
        cluster->SetYE(0.002);
      } else {
        if (_diagonal)
          cluster->SetYE(0.0008);
        else
          if (_iteration)
            cluster->SetYE(_uncertainty);
          else
            cluster->SetYE(0.001);
      }

      if (_verbose >= v_residuals) {
        std::cout << "Cluster pos " << track_pos[clusterId] << "\t";
        std::cout << cluster_mean[clusterId] << std::endl;
      }
    } // loop over clusters

    std::vector<TCluster*> clusters_clean;
    if (_diagonal) {
      // average 2 measurements into one
      for (uint pairIt = 0; pairIt < robust_clusters.size() - 1; pairIt += 2) {
        auto cluster1 = robust_clusters[pairIt+0];
        auto cluster2 = robust_clusters[pairIt+1];
        if (!cluster1 || !cluster2)
          continue;
        float x1 = geom::GetXposPad(cluster1->GetHits()[0],
                                    _invert,
                                    units::a45
                                    );
        float x2 = geom::GetXposPad(cluster2->GetHits()[0],
                                    _invert,
                                    units::a45
                                    );
        float av_x = 0.5 * (x1 + x2);

        float y1 = cluster1->GetY();
        float y2 = cluster2->GetY();

        float y1_e = cluster1->GetYE();
        float y2_e = cluster2->GetYE();

        float av_y = y1 * y2_e * y2_e + y2 * y1_e * y1_e;
        av_y /= y1_e * y1_e + y2_e * y2_e;

        auto cl = new TCluster();
        cl->SetPos(av_x, av_y, 0.);
        clusters_clean.push_back(cl);
      } // loop over pairs
    } else {
      // do nothing in case of row/column
      for (auto cluster:robust_clusters)
        clusters_clean.push_back(cluster);
    }

    if (_verbose >= v_analysis_steps) {
      std::cout << "Loop over columns done" << std::endl;
    }

    _Cols_used->Fill(Ndots);

    _sw_partial[2]->Stop();
    _sw_partial[3]->Start(false);

//******************** STEP 3 **************************************************

    TF1* fit = NULL;
    fit = _fitter->FitTrack(clusters_clean);

    if (!fit)
      continue;

    // TODO review this mess
    if (fit && _do_full_track_fit)
      fit = (TF1*)fit->Clone();

    double quality = fit->GetChisquare() / fit->GetNDF();
    _Chi2_track->Fill(quality);

    if (_verbose >= v_analysis_steps)
      std::cout << "Track fit done" << std::endl;

    TString func = fit->GetName();

    // in case of arc fitting fill the momentum histo
    // TODO review mess
    if (!_do_linear_fit && !_do_para_fit) {
      float mom = fit->GetParameter(0) * units::B * units::clight / 1.e9;
      if (func.CompareTo("circle_dn") == 0) mom *= -1.;
      _mom_reco->Fill(mom);
      _pos_reco->Fill(fit->GetParameter(2));
      _ang_reco->Fill(fit->GetParameter(1));
    }

//****************** STEP 4 ****************************************************

    TF1* fit1[Nclusters];
    for (int i = 0; i < Nclusters; ++i) {
      if (!_correction)
        fit1[i] = fit;
      else {
        fit1[i] = _fitter->FitTrack(clusters_clean, i);
        // TODO review this mess
        if (_do_full_track_fit)
          fit1[i] = (TF1*)fit1[i]->Clone();
      }
    }

    _sw_partial[3]->Stop();
    _sw_partial[4]->Start(false);

//****************** STEP 5 ****************************************************

    // second loop over columns
    for (uint clusterId = 0; clusterId < clusters_clean.size(); ++clusterId) {
      auto cluster = clusters_clean[clusterId];
      if (!cluster) continue;

      _x_av[clusterId]      = cluster->GetX();
      double track_fit_y    = fit->Eval(_x_av[clusterId]);
      double track_fit_y1   = fit1[clusterId]->Eval(_x_av[clusterId]);

      if (_verbose >= v_residuals) {
        std::cout << "Residuals id:x:cluster:track\t" << clusterId << "\t";
        std::cout << _x_av[clusterId] << "\t";
        std::cout << clusters_clean[clusterId]->GetY() << "\t";
        std::cout << track_fit_y << "\t";
        std::cout << (clusters_clean[clusterId]->GetY() - track_fit_y)*1e6;
        std::cout << std::endl;
      }

      // fill SR
      _cluster_av[clusterId] = cluster->GetY();
      _track_pos[clusterId] = track_fit_y;
      _residual[clusterId] = _cluster_av[clusterId] - track_fit_y;

      _resol_col_hist[clusterId]->Fill(_cluster_av[clusterId] - track_fit_y);
      _resol_col_hist_except[clusterId]->Fill(_cluster_av[clusterId] - track_fit_y1);
    }

// ************ STEP 6 *********************************************************

    for (uint clusterId = 0; clusterId < clusters.size(); ++clusterId) {
      a_peak_fit[clusterId] = clusters[clusterId]->GetCharge();

      // don't fill PRF for multiplicity of 1
      if (clusters[clusterId]->GetHits().size() == 1)
        continue;
      // Fill PRF
      auto robust_pads = GetRobustPadsInColumn(clusters[clusterId]->GetHits());
      int padId = 0;
      for (auto pad:robust_pads) {
        if (!pad)
          continue;

        auto q    = pad->GetQ();
        auto time = pad->GetTime();

        double x = geom::GetXposPad(pad, _invert, _diagonal ? units::a45 : 0);
        double center_pad_y = geom::GetYposPad(pad,
                                               _invert,
                                               _diagonal ? units::a45 : 0
                                               );

        double track_fit_y    = fit->Eval(x);

        /// WARNING only 9 first pads are used
        if (padId > 9)
          continue;

        _time[clusterId][padId]    = time;
        _dx[clusterId][padId]    = center_pad_y - track_fit_y;
        _qfrac[clusterId][padId] = q / a_peak_fit[clusterId];

        if (_verbose >= v_prf) {
          std::cout << "PRF fill\t" << _dx[clusterId][padId] << "\t" << _qfrac[clusterId][padId] << std::endl;
        }

        // fill PRF
        _PRF_histo->Fill( _dx[clusterId][padId],
                          _qfrac[clusterId][padId]
                          );
        _PRF_histo_col[clusterId]->Fill( _dx[clusterId][padId],
                                         _qfrac[clusterId][padId]
                                         );

        if (_multiplicity[clusterId] == 2)
          _PRF_histo_2pad->Fill( _dx[clusterId][padId],
                                 _qfrac[clusterId][padId]
                                 );
        else if (_multiplicity[clusterId] == 3)
          _PRF_histo_3pad->Fill( _dx[clusterId][padId],
                                 _qfrac[clusterId][padId]
                                 );
        else if (_multiplicity[clusterId] == 4)
          _PRF_histo_4pad->Fill( _dx[clusterId][padId],
                                 _qfrac[clusterId][padId]
                                 );

        // robust_pads are assumed sorted!!
        ++padId;
      }
    } // loop over colums
    _sw_partial[4]->Stop();
    for (int i = 0; i < Nclusters; ++i)
      if (fit1[i] && _correction) {
        delete fit1[i];
        fit1[i] = NULL;
      }

    delete fit;

    if(_test_mode) this->DrawSelectionCan(event,trackId);
    if (!_batch)
      Draw();
  } // loop over tracks

  _tree->Fill();

  if (_store_event)
    _passed_events.push_back(event->GetID());

  return true;
}

//******************************************************************************
TCanvas* SpatialResolAna::DrawSelectionCan(const TEvent* event, int trkID) {
//******************************************************************************
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

//******************************************************************************
bool SpatialResolAna::WriteOutput() {
//******************************************************************************
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
  // _PRF_graph_2pad->Fit("PRF_function_2pad", "Q", "", -0.014, 0.014);
  // _PRF_graph_3pad->Fit("PRF_function_3pad", "Q", "", fit_bound_left, fit_bound_right);
  // _PRF_graph_4pad->Fit("PRF_function_4pad", "Q", "", fit_bound_left, fit_bound_right);

  std::cout << "done" << std::endl;
  std::cout << "Process x histoes........................";

  for (auto i = 0; i < 4; ++i)
    ProfilePRF(_PRF_histo_xscan[i], _PRF_graph_xscan[i]);

  std::cout << "done" << std::endl;
  std::cout << "Process histoes..........................";

  TH1F* resol = new TH1F("resol", "", 1000, 0., 0.001);

  for (auto i = 0; i <= Nclusters - 1; ++i)
      _resol_total->Add(_resol_col_hist[i]);

  for (auto i = 1; i < Nclusters - 1; ++i) {

    TH1F* res     = _resol_col_hist[i];
    TH1F* res_e   =  _resol_col_hist_except[i];

    Double_t mean, sigma, sigma_ex;
    Double_t mean_e = 0., sigma_e = 0., sigma_ex_e = 0.;

    sigma = 0.5 * GetFWHM(res, mean);

    if (!res->Integral()) {
      if (i < 10)
        std::cout << "WARNING. SpatialResolAna::WriteOutput(). Empty residuals at " << i << std::endl;
      continue;
    }

    if (_gaussian_residuals) {
      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);
      _resol_col_hist_2pad[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad[i]->Fit("gaus", "Q");



      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cout << "WARNING. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        continue;
      }

      res_e->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);;
      _resol_col_hist_2pad_except[i]->Fit("gaus", "Q");
      _resol_col_hist_3pad_except[i]->Fit("gaus", "Q");

      TF1* func_ex = res_e->GetFunction("gaus");
      if (!func_ex) {
        std::cout << "WARNING. SpatialResolAna::WriteOutput(). Exeptional residual fit fail" << std::endl;
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
        std::cout << "WARNING. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
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
    if (!_resol_total->Integral()) {
      std::cout << "WARNING. SpatialResolAna::WriteOutput(). Empty global residual" << std::endl;
    } else {
      auto func = _resol_total->GetFunction("gaus");
      TString output = "NaN";
      if (func)
        output = TString().Itoa(1.e6*func->GetParameter(2), 10) + " um";

      std::cout << "Spatial resolution\t" << output << std::endl;
    }
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

//******************************************************************************
bool SpatialResolAna::ProfilePRF(const TH2F* PRF_h, TGraphErrors* gr) {
//******************************************************************************
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

//******************************************************************************
Double_t SpatialResolAna::GetFWHM(const TH1F* h, Double_t& mean) {
//******************************************************************************
  if (!h->Integral())
    return -1.;

  mean = h->GetMean();
  float max   = h->GetMaximum();
  float start = h->GetBinLowEdge(h->FindFirstBinAbove(max/2));
  float end   = h->GetBinLowEdge(h->FindLastBinAbove(max/2)) +
  h->GetBinWidth(h->FindLastBinAbove(max/2));

  return end - start;
}

//******************************************************************************
TF1* SpatialResolAna::InitializePRF(const TString name) {
//******************************************************************************
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

  double co = 0.85436010;
  double a2 = -1919.1462;
  double a4 = 17575063.;
  double b2 = 627.88307;
  double b4 = 1.0339875e+09;

  func->SetParameters(co, a2, a4, b2, b4);

  return func;
}

bool SpatialResolAna::Draw() {
  std::cout << "Draw event " << _ev << std::endl;
  TCanvas c("c_SR", "c_SR", 800, 600);
  TGraphErrors* gr = new TGraphErrors();
  TGraphErrors* gr_f = new TGraphErrors();
  TGraphErrors* gr_c = new TGraphErrors();
  auto scale = 1.e3;
  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_cluster_av[colIt] == -999)
      continue;

    gr->SetPoint(gr->GetN(), scale*_x_av[colIt], scale*_cluster_av[colIt]);
    gr_f->SetPoint(gr_f->GetN(), scale*_x_av[colIt], scale*_track_pos[colIt]);
  }

  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_clust_pos[colIt] == -999)
      continue;

    gr_c->SetPoint(gr_c->GetN(), scale*_x[colIt], scale*_clust_pos[colIt]);
  }

  gr_c->SetTitle("Event " + TString().Itoa(_ev, 10));
  gr_c->GetYaxis()->SetTitle("Reference Y, [mm]");
  gr_c->GetXaxis()->SetTitle("Reference X, [mm]");
  gPad->SetGrid();
  gr_c->SetMarkerStyle(kPlus);
  gr_c->Draw("ap");
  gr->Draw("same p");


  gr_f->SetLineColor(kRed);
  gr_f->Draw("same l");
  c.Update();
  c.WaitPrimitive();
  return true;
}

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna(argc, argv);
  if (!ana->Initialize())               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}

