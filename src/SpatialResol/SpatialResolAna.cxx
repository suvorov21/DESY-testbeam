/** @cond */
#include "TVector3.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
/** @endcond */

#include "SpatialResolAna.hxx"
#include "line.hxx"

//******************************************************************************
SpatialResolAna::SpatialResolAna():
  AnalysisBase(),
  _prev_iter_file(nullptr),
  _prev_iter_name(TString("")),
  _iteration(0),
  _correction(false),
  _tree(nullptr),
  _do_full_track_fit(false),
  _do_separate_pad_fit(false),
  _gaussian_residuals(true),
  _charge_uncertainty(true),
  _prf_scale_axis(nullptr),
  _x_scan_axis(nullptr) {
//******************************************************************************
}

//******************************************************************************
bool SpatialResolAna::ReadCLI(int argc, char **argv) {
//******************************************************************************
  if (!AnalysisBase::ReadCLI(argc, argv))
    return false;
  char * pEnd;

  for (auto i = 1; i < argc; ++i) {
    std::string argString = argv[i];
    _prev_iter_name = lookForOption(argString, i == argc - 1 ? "" : argv[i + 1],
                                    {"--prev"}, _prev_iter_name.Data());
    _iteration = (int)strtol(lookForOption(
        argString, i == argc - 1 ? "" : argv[i + 1], {"-t", "--iter"}, std::to_string(_iteration)).c_str(), &pEnd, 10);

    std::string corrOpt = lookForOption(
        argString, "t", {"-c", "--corr"}, "f");

    if (corrOpt == "t")
      _correction = true;
  }

  return true;
}

//******************************************************************************
bool SpatialResolAna::Initialize(int argc, char** argv) {
//******************************************************************************
  std::cout << "*****************************************" << std::endl;
  std::cout << "***   Spatial resolution analysis    ****" << std::endl;
  std::cout << "*****************************************" << std::endl;

  AnalysisBase::Initialize(argc, argv);

  std::cout << "Batch mode    :   " << _batch         << std::endl;
  std::cout << "Verbosity     :   " << _verbose       << std::endl;
  std::cout << "Iteration     :   " << _iteration     << std::endl;
  std::cout << "Debug         :   " << _test_mode     << std::endl;

  std::cout << "Initializing spatial resolution ana......";

  gErrorIgnoreLevel = kSysError;

  _uncertainty_vs_prf_gr_prev = nullptr;
  _uncertainty_vs_prf_histo   = nullptr;

  _prf_function = InitializePRF("PRF_function", _prf_free_centre);
  // _prf_function_2pad = InitializePRF("PRF_function_2pad");
  // _prf_function_3pad = InitializePRF("PRF_function_3pad");
  // _prf_function_4pad = InitializePRF("PRF_function_4pad");

  if (_do_full_track_fit)
    _prf_function->FixParameter(0, 1.);

  // Initialize graph for PRF profiling
  _prf_graph = new TGraphErrors();
  _prf_graph->SetName("PRF_graph");
  _prf_graph_2pad = new TGraphErrors();
  _prf_graph_2pad->SetName("PRF_graph_2pad");
  _prf_graph_3pad = new TGraphErrors();
  _prf_graph_3pad->SetName("PRF_graph_3pad");
  _prf_graph_4pad = new TGraphErrors();
  _prf_graph_4pad->SetName("PRF_graph_4pad");

  // load information from previous iteration
  if (_iteration) {
    if (_prev_iter_name.Length() == 0) {
      _prev_iter_name = _file_out_name;
      _prev_iter_name   = _prev_iter_name(0, _prev_iter_name.Index("iter"));
      _prev_iter_name  += "iter";
      _prev_iter_name  += TString::Itoa(_iteration - 1, 10);
      _prev_iter_name  += ".root";
    }

    _prev_iter_file = new TFile(_prev_iter_name.Data(), "READ");
    if (!_prev_iter_file->IsOpen()) {
      std::cerr << "ERROR! " << __func__ << std::endl;
      std::cerr << "File from previous iteration is not found" << std::endl;
      std::cerr << "File name: " << _prev_iter_name << std::endl;
      exit(1);
    }

    // READ PRF
    TH2F* histo_prev = (TH2F*)_prev_iter_file->Get("PRF_histo");
    auto tree = (TTree*)_prev_iter_file->Get("outtree");
    if (!histo_prev) {
      histo_prev = new TH2F("PRF_histo_tmp2","", prf_bin, prf_min, prf_max, 150,0.,1.5);
      tree->Project("PRF_histo_tmp2", "qfrac:dx");
    }
    histo_prev->SetName("prev_hsto");
    if (!ProfilePRF(histo_prev, _prf_graph)) {
      std::cerr << "ERROR! " << __func__  << std::endl;
      std::cerr << "PRF can not be profiled" << std::endl;
      exit(1);
    }
    // kind of magic that works.
    // More fits better result
    // seems that the tolerance is wronng somewhere in ROOT
    for (auto i = 0; i < 3; ++i)
      _prf_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);
    // TODO think about memory control
    _prf_function = (TF1*)_prf_graph->GetFunction("PRF_function")->Clone("PRF_function");
    auto uncertainty_graph = (TH1F*)_prev_iter_file->Get("resol_total");

    if (!_prf_function || !uncertainty_graph) {
      std::cerr << "ERROR. " << __func__ ;
      std::cerr << "PRF function or resolution is not specified" << std::endl;
      std::cerr << "Search in " << _prev_iter_name << std::endl;
      exit(1);
    }

    if (_clustering->n_pads > 0 && _individual_column_PRF) {
      std::cerr << "ERROR. Conflicting options" << std::endl;
      exit(1);
    }

    // Read PRF for complicated patterns
    if (_clustering->n_pads > 1) {
      _prf_function_arr = new TF1*[3];
      for (auto rest = 0; rest < _clustering->n_pads; ++rest) {
        _prf_function_arr[rest] = InitializePRF("PRF_function_tmp", true);

        TH2F* tmp = new TH2F("PRF_histo_tmp","", prf_bin, prf_min, prf_max, 150,0.,1.5);
        TString s = TString().Itoa(_clustering->n_pads, 10);
        TString r = TString().Itoa(rest, 10);
        tree->Project("PRF_histo_tmp", "qfrac:dx", "abs(pad_x%" + s + ") == " + r);
        auto gr_tmp = new TGraphErrors();
        ProfilePRF(tmp, gr_tmp);
        for (auto i = 0; i < 3; ++i)
          gr_tmp->Fit("PRF_function_tmp", "Q", "", fit_bound_left, fit_bound_right);

        _prf_function_arr[rest] = (TF1*)gr_tmp->GetFunction("PRF_function_tmp")->Clone(Form("PRF_function_arr_%i", rest));
      }
    }

    if (_individual_column_PRF) {

      _prf_function_arr = new TF1*[36];
      for (auto colId = 0; colId < geom::GetNColumn(_invert); ++colId) {
        _prf_function_arr[colId] = InitializePRF("PRF_function_tmp", _prf_free_centre);

        TH2F* tmp = new TH2F("PRF_histo_tmp","", prf_bin, prf_min, prf_max, 150,0.,1.5);
        TString s = TString().Itoa(_clustering->n_pads, 10);
        TString r = TString().Itoa(colId, 10);
        tree->Project("PRF_histo_tmp", "qfrac:dx", "pad_x == " + r);
        auto gr_tmp = new TGraphErrors();
        ProfilePRF(tmp, gr_tmp);
        for (auto i = 0; i < 3; ++i)
          gr_tmp->Fit("PRF_function_tmp", "Q", "", fit_bound_left, fit_bound_right);

        _prf_function_arr[colId] = (TF1*)gr_tmp->GetFunction("PRF_function_tmp")->Clone(Form("PRF_function_arr_%i", colId));
      }
    }

    // Read PRF in time
    auto histo_prev_t = (TH2F*)_prev_iter_file->Get("PRF_histo_time");
    histo_prev_t->SetName("prev_hsto_time");
    auto Nbins = histo_prev_t->GetYaxis()->GetNbins();
    _prf_time_e = new TH1F("prf_e", "",
                           Nbins,
                           histo_prev_t->GetYaxis()->GetBinLowEdge(1),
                           histo_prev_t->GetYaxis()->GetBinLowEdge(Nbins) + histo_prev_t->GetYaxis()->GetBinWidth(Nbins)
                           );

    _prf_time_error = new TGraphErrors();
    if (!_prf_time_error || !ProfilePRF_X(histo_prev_t, _prf_time_error, _prf_time_e)) {
      std::cerr << "ERROR. SpatialResolAna::Initialize().";
      std::cout << "PRF time function is not specified" << std::endl;
      std::cerr << "Search in " << _prev_iter_name << std::endl;
      exit(1);
    }
    histo_prev_t->Fit("pol2", "Q");
    _prf_time_error->Fit("pol2", "Q", "", -0.015, -0.005);
    _prf_time_error->Fit("pol2", "Q", "", -0.015, -0.005);
    _prf_time_func = _prf_time_error->GetFunction("pol2");

    Double_t mean, sigma;
    sigma = 0.5 * GetFWHM(uncertainty_graph, mean);
    uncertainty_graph->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);
    _uncertainty = (Float_t)uncertainty_graph->GetFunction("gaus")->GetParameter(2);

    TGraphErrors* temp = nullptr;
    if (_do_separate_pad_fit)
      temp = (TGraphErrors*)_prev_iter_file->Get("uncertainty_vs_prf_gr");
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

    // read event list passed through reconstruction+selection at previous iteration
    Int_t read_var;
    auto event_tree = (TTree*)_prev_iter_file->Get("EventTree");
    event_tree->SetBranchAddress("PassedEvents",    &read_var);
    std::vector<Int_t> vec;
    vec.clear();
    for (auto i = 0; i < event_tree->GetEntries(); ++i) {
      event_tree->GetEntry(i);
      vec.push_back(read_var);
    }
    this->SetEventList(vec);

  } // if iteration

  _qulity_ratio   = new TH1F("quality_ratio",
    "Quality ratio arc/linear", 100, 0., 2.);
  _chi2_ratio     = new TH1F("chi_ratio",
    "Chi2 ratio arc/linear", 100, 0., 2.);

  _mom_reco   = new TH1F("mom_reco", "", 2000, -6., 6.);
  _pos_reco   = new TH1F("pos_reco", "", 8000, -0.2, 0.2);
  _ang_reco   = new TH1F("ang_reco", "", 3000, -0.3, 0.3);

  // Initialise histoes and graphs
  _prf_histo = new TH2F("PRF_histo","", prf_bin, prf_min, prf_max, 150,0.,1.5);

  // PRF for different multiplicities
  auto dir_prf_mult = _file_out->mkdir("prf_mult");
  _output_vector.push_back(dir_prf_mult);

  // PRF in time
  _prf_time = new TH2F("PRF_histo_time","", prf_bin, prf_min, prf_max, 100, -20., 80.);
  _output_vector.push_back(_prf_time);

  _file_out->cd();
  _tree = new TTree("outtree", "");
  _tree->Branch("ev",           &_ev);
  _tree->Branch("dEdx",         &_dEdx);

  _tree->Branch("angle_yz",     &_angle_yz);
  _tree->Branch("angle_xy",     &_angle_xy);

  _tree->Branch("rob_clusters", &_rob_clusters);

  _tree->Branch("quality",      &_quality);
  _tree->Branch("mom",          &_mom);
  _tree->Branch("sina",         &_sin_alpha);
  _tree->Branch("offset",       &_offset);

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
  _tree->Branch("residual_corr",
                &_residual_corr,
                TString::Format("residual_corr[%i]/F", Nclusters)
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
  _tree->Branch("pad_charge",
                &_pad_charge,
                TString::Format("_pad_charge[%i][10]/I", Nclusters)
                );
  _tree->Branch("pad_time",
                &_pad_time,
                TString::Format("_pad_time[%i][10]/I", Nclusters)
                );

  _tree->Branch("wf_width",
                &_wf_width,
                TString::Format("_wf_width[%i][10]/I", Nclusters)
                );

  _tree->Branch("wf_fwhm",
                &_wf_fwhm,
                TString::Format("_wf_fwhm[%i][10]/I", Nclusters)
                );

  _tree->Branch("pad_x",
                &_pad_x,
                TString::Format("_pad_x[%i][10]/I", Nclusters)
                );
  _tree->Branch("pad_y",
                &_pad_y,
                TString::Format("_pad_y[%i][10]/I", Nclusters)
                );

  if (_to_store_wf) _tree->Branch("pad_wf_q",
                &_pad_wf_q,
                TString::Format("_pad_wf_q[%i][10][520]/I", Nclusters)
                );

  // WARNING TEMP
  _tree->Branch("fit_up",
                &_fit_up,
                TString::Format("fit_up[%i]/F", Nclusters)
                );

  _tree->Branch("fit_bt",
                &_fit_bt,
                TString::Format("fit_bt[%i]/F", Nclusters)
                );

  _output_vector.push_back(_tree);

  _hdEdx  = new TH1F("dEdx","",300,0,5000);
  _output_vector.push_back(_hdEdx);

  _prf_histo_2pad = new TH2F("PRF_histo_2pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _prf_histo_3pad = new TH2F("PRF_histo_3pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);
  _prf_histo_4pad = new TH2F("PRF_histo_4pad",
    "", prf_bin, prf_min, prf_max, 150,0.,1.5);

  dir_prf_mult->Append(_prf_histo_2pad);
  dir_prf_mult->Append(_prf_histo_3pad);
  dir_prf_mult->Append(_prf_histo_4pad);

  // dir_prf_mult->Append(_prf_function_2pad);
  // dir_prf_mult->Append(_prf_function_3pad);
  // dir_prf_mult->Append(_prf_function_4pad);

  dir_prf_mult->Append(_prf_graph_2pad);
  dir_prf_mult->Append(_prf_graph_3pad);
  dir_prf_mult->Append(_prf_graph_4pad);

  for (auto i = 0; i < Nclusters; ++i)
    _prf_histo_col[i] = new TH2F(Form("PRF_histo_col_%i", i),
      "", prf_bin, prf_min, prf_max, 150,0.,1.5);

  for (auto i = 0; i < 4; ++i) {
    _prf_histo_xscan[i] = new TH2F(Form("PRF_histo_pad_%i", i),
      "", prf_bin, prf_min, prf_max, 150,0.,1.5);
    _prf_graph_xscan[i] = new TGraphErrors();
    _prf_graph_xscan[i]->SetName(Form("PRF_graph_pad_%i", i));
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

  _chi2_track = new TH1F("Chi2_Track", "", 1000, 0., 100.);
  _cols_used  = new TH1F("cols_used", "", 40, 0., 40.);

  // schedule the output for writing
  _output_vector.push_back(_prf_function);
  _output_vector.push_back(_prf_histo);
  _output_vector.push_back(_prf_graph);

  _output_vector.push_back(_residual_sigma);
  _output_vector.push_back(_residual_mean);

  _output_vector.push_back(_residual_sigma_biased);
  _output_vector.push_back(_residual_sigma_unbiased);

  _output_vector.push_back(_chi2_track);
  _output_vector.push_back(_cols_used);
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
        _fit_quality_plots[id][j] = new TH1F(Form("pad_fit_q_%i_%i", j, id),
         "", 300, 0., 30.);
        _output_vector.push_back(_fit_quality_plots[id][j]);
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
  for (auto & j : _prf_histo_col) {
    dir_prf->Append(j);
  }

  _passed_events.clear();

  std::cout << "done" << std::endl;
  if (_verbose >= v_event_number) {
    std::cout << "\t PRF(x) = " << _prf_function->GetFormula()->GetExpFormula();
    std::cout << "  with ";
    for (auto i = 0; i < _prf_function->GetNpar(); ++i)
      std::cout << "  " << _prf_function->GetParameter(i) << ",";
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
                             _verbose,
                             _iteration,
                             _prf_function,
                             _prf_graph,
                             fit_bound_right,
                             _charge_uncertainty,
                             _prf_time_func,
                             _prf_time_e,
                             _clustering->angle
                             );

  if (_clustering->n_pads > 1 && _iteration) {
    std::cout << _prf_function_arr[0]->GetName() << std::endl;
    std::cout << _prf_function_arr[1]->GetName() << std::endl;
    _fitter->SetPRFarr(_prf_function_arr, _clustering->n_pads);
    _fitter->SetComplicatedPatternPRF(true);
  }

  if (_individual_column_PRF && _iteration) {
    _fitter->SetPRFarr(_prf_function_arr, 36);
    _fitter->SetIndividualPRF(true);
  }

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

  auto track_hits = event->GetUsedHits();
  if (track_hits.empty())
    return false;

  _ev = event->GetID();

  _quality    = -999.;
  _mom        = -999.;
  _sin_alpha  = -999.;
  _offset     = -999.;
  _rob_clusters = -999;

  int cluster_N[Nclusters];
  double track_pos[Nclusters];
  double cluster_mean[Nclusters];
  float charge_max[Nclusters];
  double a_peak_fit[Nclusters];

  int Ndots = 0;

  // At the moment ommit first and last column
  // first loop over columns
  _sw_partial[2]->Start(false);

  // reset tree values
  for (auto colId = 0; colId < Nclusters; ++colId) {
    _multiplicity[colId]  = -999;
    _charge[colId]        = -999;
    _residual[colId]      = -999;
    _residual_corr[colId] = -999;
    _clust_pos[colId]     = -999;
    _track_pos[colId]     = -999;
    _x[colId]             = -999;
    _x_av[colId]          = -999;
    // _cluster_av[colId]    = -999;
    _dEdx               = -999;

    // WARNING TMP
    _fit_up[colId]        = -999.;
    _fit_bt[colId]        = -999.;

    for (auto padId = 0; padId < 10; ++padId) {
      _dx[colId][padId]    = -999;
      _time[colId][padId]  = -999;
      _qfrac[colId][padId] = -999;

      //dEdx part
      _pad_charge[colId][padId] = -999;
      _pad_time[colId][padId] = -999;
      _pad_x[colId][padId] = -999;
      _pad_y[colId][padId] = -999;

      _wf_width[colId][padId] = -999;
      _wf_fwhm[colId][padId] = -999;

     for (auto nSmp = 0; nSmp < 520; ++nSmp) {
    //_pad_wf_t[colId][padId][nSmp] = -999;
    _pad_wf_q[colId][padId][nSmp] = -999;

    	}

    }

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
  clusters = ClusterTrack(track_hits);

  if (_verbose >= v_analysis_steps)
    std::cout << "Clusterization done " << clusters.size() <<  std::endl;

  // selection
  if (!sel::CrossingTrackSelection(clusters,
                                   _max_mult,
                                   _cut_gap,
                                   _max_phi,
                                   _max_theta,
                                   _invert,
                                   _verbose
                                   ))
    return false;

  std::vector<double> fit_v = sel::GetFitParams(clusters, _invert);
  std::vector<double> fit_xz = sel::GetFitParamsXZ(clusters, _invert);

  _angle_xy = (Float_t)fit_v[2];
  _angle_yz = (Float_t)fit_xz[2] * sel::v_drift_est;

  // if not a column clustering
  if (_clustering->n_pads > 0) {
    if (clusters.size() < 5)
      return false;
    // clean first and last cluster
    sort(clusters.begin(), clusters.end(),
         [&](TCluster* cl1, TCluster* cl2){
            return cl1->GetX() < cl2->GetX();
          });
    clusters.erase(clusters.begin(), clusters.begin() + 1);
    clusters.erase(clusters.end()-1);
  }
  // cut on the number of clusters
  if (clusters.size() < uint(_min_clusters))
    return false;
  // truncation
  auto robust_clusters = GetRobustClusters(clusters);

  if (_verbose >= v_analysis_steps)
    std::cout << "clearing done, columns\t" << robust_clusters.size() << std::endl;

  _rob_clusters = (int)robust_clusters.size();

  /// decide that track is accepted by selection
  _store_event = true;
// *******************  STEP 2 *************************************************

  if (_verbose >= v_analysis_steps)
    std::cout << "start cluster fit" << std::endl;

  std::vector <double> QsegmentS; QsegmentS.clear();
  std::vector <int> pad_wf_v; //WF

  //if (robust_clusters.size() < 30) continue;
        //continue;

  for (uint clusterId = 0; clusterId < robust_clusters.size(); ++clusterId) {
    auto cluster = robust_clusters[clusterId];
    if (!(&cluster[0]))
      continue;

    // loop over rows
    auto robust_pads = GetRobustPadsInCluster(cluster->GetHits());
    _multiplicity[clusterId] = (int)robust_pads.size();

    _clust_pos[clusterId] = 0.;
    _charge[clusterId] = 0;

    int colQ = 0;
    int padId = 0;

    auto pad_id = -1;
    for (auto pad:robust_pads) {
      ++pad_id;
      if (!pad)
        continue;

      // treat cross-talk
      // not for the first pad
      if (pad != robust_pads[0] && _cross_talk_treat != def) {
        auto dt = abs(pad->GetTime() - robust_pads[0]->GetTime());
        auto qfrac = 1. * pad->GetQ() / robust_pads[0]->GetQ();
        // cross talk selection
        if (_cross_talk_treat == suppress && dt < 4 && qfrac < 0.08) {
          pad->SetQ(0);
          for (auto time = pad->GetTime() + 4; time < 510; ++time) {
            if (pad->GetADC(time) > pad->GetQ()) {
              pad->SetTime(time);
              pad->SetQ(pad->GetADC(time));
            }
          } // over time
          if (pad->GetQ() == 0) {
            robust_pads.erase(
                              robust_pads.begin() + pad_id,
                              robust_pads.begin() + pad_id + 1
                              );
            --pad_id;
          }
        } // cross-talk suppression

        if (_cross_talk_treat == cherry_pick) {
          pad->SetQ(0);
          for (
               auto time = robust_pads[0]->GetTime() - 4;
               time < robust_pads[0]->GetTime() + 4;
               ++time
               ) {
            if (pad->GetADC(time) > pad->GetQ() &&
                pad->GetADC(time) > pad->GetADC(robust_pads[0]->GetTime() + 5)) {
              pad->SetTime(time);
              pad->SetQ(pad->GetADC(time));
            }
          } // over time
          if (pad->GetQ() == 0) {
            robust_pads.erase(
                              robust_pads.begin() + pad_id,
                              robust_pads.begin() + pad_id + 1
                              );
            --pad_id;
          }
        } // cross-talk cherry picking
      } // cross-talk block

      _clust_pos[clusterId] += (float) pad->GetQ() * geom::GetYposPad(pad,
                                                                      _invert,
                                                                      _clustering->angle);
      _charge[clusterId] += pad->GetQ();


      //dEdx part
    	pad_wf_v.clear();

      colQ+= pad->GetQ();

      _pad_charge[clusterId][padId] = pad->GetQ();
      _pad_time[clusterId][padId] = pad->GetTime();
      _pad_x[clusterId][padId] = pad->GetCol(_invert);
      _pad_y[clusterId][padId] = pad->GetRow(_invert);
      _wf_width[clusterId][padId] = pad->GetWidth();
      _wf_fwhm[clusterId][padId] = pad->GetFWHM();

      //if (_pad_charge[clusterId][padId] <= 0) continue;

       if (_to_store_wf) {
           for (int tz = 0; tz < geom::Nsamples; ++tz) {
               _pad_wf_q[clusterId][padId][tz] = pad->GetADC(tz);
           }
    	}
      ++padId;
    } // loop over pads
    _clust_pos[clusterId] /= (Float_t)_charge[clusterId];
    _x[clusterId] = cluster->GetX();

    double CoC =  _clust_pos[clusterId];

    ++Ndots;

    if (_multiplicity[clusterId] > 1 && _iteration > 0) {
      track_pos[clusterId] = _fitter->FitCluster(robust_pads,
                                                 _charge[clusterId],
                                                 _clust_pos[clusterId]
                                                 );
      _clust_pos[clusterId] = (float_t)track_pos[clusterId];
      // WARNING temp
      auto it_main = std::max_element((*cluster).begin(), (*cluster).end(),
                                      [](const THit* n1, const THit* n2) { return n1->GetQ() < n2->GetQ(); });
      auto main_row = (*it_main)->GetRow(_invert);
      auto it_up = std::find_if((*cluster).begin(), (*cluster).end(),
                                    [&](const THit* h1) { return h1->GetRow(_invert) == main_row + 1; });
      auto it_bt = std::find_if((*cluster).begin(), (*cluster).end(),
                                    [&](const THit* h1) { return h1->GetRow(_invert) == main_row - 1; });

      // for (auto ii = 0; ii < cluster->GetSize(); ++ii)
      //   std::cout << (*cluster)[ii]->GetRow(_invert) << "   ";
      // std::cout << std::endl;

      if (it_up != (*cluster).end() && it_bt != (*cluster).end()) {
        // std::cout << main_row << "\t" << (*it_up)->GetRow() << "\t" << (*it_bt)->GetRow() << std::endl;
        Float_t r_up = (Float_t)(*it_up)->GetQ() / (Float_t)(*it_main)->GetQ();
        Float_t r_bt = (Float_t )(*it_bt)->GetQ() / (Float_t)(*it_main)->GetQ();
        if (r_up > _prf_function->GetMinimum() && r_bt > _prf_function->GetMinimum()) {
          // r_up = PRF(x-x_pad) / PRF(x - x_pad+1)
          auto minimisator = [&](const Double_t *par) {
            if (abs(par[1]) < 1e-6)
              return abs(r_up - _prf_function->Eval(geom::GetYposPad((*it_up)) - par[0]) / _prf_function->Eval(geom::GetYposPad((*it_main)) - par[0]));
            else
              return abs(r_bt - _prf_function->Eval(geom::GetYposPad((*it_bt)) - par[0]) / _prf_function->Eval(geom::GetYposPad((*it_main)) - par[0]));
          };
          ROOT::Math::Functor fcn_cluster(minimisator,2);
          ROOT::Fit::Fitter  fitter_cluster;

          double pStart[2] = {geom::GetYposPad((*it_main)), true};
          fitter_cluster.SetFCN(fcn_cluster, pStart);
          bool ok = fitter_cluster.FitFCN();
          (void)ok;
          const ROOT::Fit::FitResult & result_cluster = fitter_cluster.Result();
          if (ok && r_up > 0.04)
            _fit_up[clusterId] = (Float_t)result_cluster.GetParams()[0];

          // if (_fit_up[clusterId] > -0.009 && _fit_up[clusterId] < -0.0087)
          //   std::cout << r_up << "\t" << _fit_up[clusterId] << "\t" << _prf_function->Eval(geom::GetYposPad((*it_up)) - _fit_up[clusterId]) / _prf_function->Eval(geom::GetYposPad((*it_main)) - _fit_up[clusterId]) << std::endl;

          pStart[0] = geom::GetYposPad((*it_main));
          pStart[1] = false;
          fitter_cluster.SetFCN(fcn_cluster, pStart);
          ok = fitter_cluster.FitFCN();
          const ROOT::Fit::FitResult & result_cluster2 = fitter_cluster.Result();
          if (ok && r_bt > 0.04)
            _fit_bt[clusterId] = result_cluster2.GetParams()[0];
        }
      }
    } else {
      track_pos[clusterId] = _clust_pos[clusterId];
    }

    cluster->SetY(track_pos[clusterId]);
    cluster->SetCharge(_charge[clusterId]);

    if (_verbose >= v_fit_details) {
      for (auto pad:*cluster) {
        std::cout << pad->GetRow(_invert) <<  " : " << pad->GetCol(_invert) << "\t";
      }
      std::cout << std::endl;
    }

    if (_verbose >= v_fit_details)
      std::cout << "X:CoC:Fit\t" << cluster->GetX() << "\t" << CoC << "\t" << cluster->GetY() << std::endl;


    // TODO review this mess
    if (cluster->GetSize() == 1) {
      cluster->SetYE(0.002);
    } else {
      // if not a column clustering
      if (_clustering->n_pads > 0)
        cluster->SetYE(0.0008);
      else {
        /** DEV VERSION */
        // auto half_size = 0.5 * geom::dy;
        // auto dx_from_side = fmod(abs(cluster->GetY()), geom::dy);
        // if (dx_from_side > half_size)
        //   dx_from_side -= half_size;
        // else
        //   dx_from_side = half_size - dx_from_side;
        // auto dy = 200 - 100 * dx_from_side / half_size;
        // dy *= 1e-6;
        // cluster->SetYE(dy);
        /** END OF DEV */

        /** OLD VERSION below */
        if (_iteration)
          cluster->SetYE(_uncertainty);
        else
          cluster->SetYE(0.001);

      }
    }

    if (_verbose >= v_residuals) {
      std::cout << "Cluster pos " << track_pos[clusterId] << "\t";
      std::cout << cluster_mean[clusterId] << std::endl;
    }

      if (colQ) QsegmentS.push_back(colQ);

  } // loop over clusters

  std::vector<TCluster*> clusters_clean;

  for (auto cluster:robust_clusters)
    clusters_clean.push_back(cluster);

  if (_verbose >= v_analysis_steps) {
    std::cout << "Loop over columns done" << std::endl;
  }

    //dEdx calculations
    double alpha = 0.7;
    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    float CT= totQ / (alpha * QsegmentS.size());

    //_npoints = QsegmentS.size();
    _dEdx = CT;
    _hdEdx->Fill(CT);

  _cols_used->Fill(Ndots);

  _sw_partial[2]->Stop();
  _sw_partial[3]->Start(false);

//******************** STEP 3 **************************************************

  TF1* fit = NULL;
  fit = _fitter->FitTrack(clusters_clean);
  _track_fit_func = fit;

  if (!fit)
    return false;

  // TODO review this mess
  if (fit && _do_full_track_fit)
    fit = (TF1*)fit->Clone();

  _quality = fit->GetChisquare() / fit->GetNDF();
  _chi2_track->Fill(_quality);

  if (_verbose >= v_analysis_steps)
    std::cout << "Track fit done" << std::endl;

  TString func = fit->GetName();

  // in case of arc fitting fill the momentum histo
  // TODO review mess
  if (!_do_linear_fit && !_do_para_fit) {
    float mom = 1./fit->GetParameter(0) * units::B * units::clight / 1.e9;
    // if (func.CompareTo("circle_dn") == 0) mom *= -1.;
    _mom_reco->Fill(mom);
    _pos_reco->Fill(fit->GetParameter(2));
    _ang_reco->Fill(fit->GetParameter(1));

    _mom        = mom;
    _sin_alpha  = fit->GetParameter(1);
    _offset     = fit->GetParameter(2);
  }

  if (_do_para_fit) {
    // pol2 := [p0]+[p1]*x+[p2]*pow(x,2)

    double x1 = clusters_clean[0]->GetX();
    TVector3 start(x1, fit->Eval(x1), 0.);
    double x2 = (*(clusters_clean.end()-1))->GetX();
    TVector3 end(x2, fit->Eval(x2), 0.);
    TLine_att line(start, end);

    double a = fit->GetParameter(2);;
    double b = fit->GetParameter(1);
    // double c = fit->GetParameter(0);

    double x0 = - b * (x2-x1) + a*x2*x2 + b*x2 - a*x1*x1 - b*x1 ;
    x0 /= 2*a*(x2-x1);
    TVector3 max_sag(x0, fit->Eval(x0), 0.);
    double sag = line.GetDistVec(max_sag).Mag();
    double L = (end - start).Mag();

    if (abs(sag) > 1e-10) {
      _mom = L*L /8/sag + sag / 2;
      _mom *= units::B * units::clight / 1.e9;
      if (a < 0)
        _mom *= -1;
    }

    _sin_alpha = TMath::Sin(TMath::ATan(fit->Derivative(x1)));
    _offset    = start.Y();

    if (_verbose > v_analysis_steps) {
      std::cout << "start:end:max\t" << start.X() << ", " << start.Y();
      std::cout << "\t" << end.X() << ", " << end.Y();
      std::cout << "\t" << max_sag.X() << ", " << max_sag.Y() << std::endl;
      std::cout << "Line at x0\t" << line.EvalX(max_sag.X()).Y() << std::endl;
      std::cout << "Sagitta:\t" << sag << "\tLength:\t" << L << std::endl;
      std::cout << "mom:\t" << _mom << std::endl;
    }
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
      std::cout << cluster->GetY() << "\t";
      std::cout << track_fit_y << "\t";
      std::cout << (cluster->GetY() - track_fit_y)*1e6;
      std::cout << std::endl;
    }

    // fill SR
    _clust_pos[clusterId] = cluster->GetY();
    _track_pos[clusterId] = track_fit_y;
    _residual[clusterId] = _clust_pos[clusterId] - track_fit_y;
    _residual_corr[clusterId] = _clust_pos[clusterId] - track_fit_y1;

    _resol_col_hist[clusterId]->Fill(_clust_pos[clusterId] - track_fit_y);
    _resol_col_hist_except[clusterId]->Fill(_clust_pos[clusterId] - track_fit_y1);
  }

// ************ STEP 6 *********************************************************

  for (uint clusterId = 0; clusterId < clusters.size(); ++clusterId) {
    auto cluster = clusters[clusterId];
    a_peak_fit[clusterId] = cluster->GetCharge();

    // don't fill PRF for multiplicity of 1
    if (cluster->GetSize() == 1)
      continue;
    // Fill PRF
    auto robust_pads = GetRobustPadsInCluster(cluster->GetHits());
    int padId = 0;

    for (auto pad:robust_pads) {
      if (!pad)
        continue;

      auto q    = pad->GetQ();
      auto time = pad->GetTime();

      double x = geom::GetXposPad(pad, _invert, _clustering->angle);
      // double x = cluster->GetX();
      double center_pad_y = geom::GetYposPad(pad,
                                             _invert,
                                             _clustering->angle
                                             );

      double track_fit_y    = fit->Eval(x);

      /// WARNING only 10 first pads are used
      if (padId > 9)
        continue;

      // Pad characteristics
      // time, position wrt track, charge, col/row
      _time[clusterId][padId]  = time;
      _dx[clusterId][padId]    = center_pad_y - track_fit_y;
      _qfrac[clusterId][padId] = q / a_peak_fit[clusterId];

      if (_verbose >= v_prf) {
        std::cout << "PRF fill\t" << _dx[clusterId][padId];
        std::cout << "\t" << _qfrac[clusterId][padId];
        std::cout << "\t" << _time[clusterId][padId] - _time[clusterId][0] << std::endl;
      }

      // fill PRF
      _prf_histo->Fill( _dx[clusterId][padId],
                        _qfrac[clusterId][padId]
                        );
      _prf_histo_col[clusterId]->Fill( _dx[clusterId][padId],
                                       _qfrac[clusterId][padId]
                                       );

      if (_multiplicity[clusterId] == 2)
        _prf_histo_2pad->Fill( _dx[clusterId][padId],
                               _qfrac[clusterId][padId]
                               );
      else if (_multiplicity[clusterId] == 3)
        _prf_histo_3pad->Fill( _dx[clusterId][padId],
                               _qfrac[clusterId][padId]
                               );
      else if (_multiplicity[clusterId] == 4)
        _prf_histo_4pad->Fill( _dx[clusterId][padId],
                               _qfrac[clusterId][padId]
                               );

      if (padId > 0)
        _prf_time->Fill(_dx[clusterId][padId], _time[clusterId][padId] - _time[clusterId][0]);

      // robust_pads are assumed sorted!!
      ++padId;
    }
  } // loop over columns
  _sw_partial[4]->Stop();
  for (int i = 0; i < Nclusters; ++i)
    if (fit1[i] && _correction) {
      delete fit1[i];
      fit1[i] = NULL;
    }

  if(_test_mode) this->DrawSelectionCan(event);
  if (!_batch)
    Draw();
  delete fit;

  _tree->Fill();

  if (_store_event)
    _passed_events.push_back(event->GetID());

  return true;
}

//******************************************************************************
bool SpatialResolAna::WriteOutput() {
//******************************************************************************
  if (!_file_out)
    return true;

  std::cout << "PRF profiling............................";

  ProfilePRF(_prf_histo, _prf_graph);
  ProfilePRF(_prf_histo_2pad, _prf_graph_2pad);
  ProfilePRF(_prf_histo_3pad, _prf_graph_3pad);
  ProfilePRF(_prf_histo_4pad, _prf_graph_4pad);

  std::cout << "done" << std::endl;
  std::cout << "PRF fit..................................";

  // MAGIC
  for (auto i = 0; i < 3; ++i)
    _prf_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);
  // _prf_graph_2pad->Fit("PRF_function_2pad", "Q", "", -0.014, 0.014);
  // _prf_graph_3pad->Fit("PRF_function_3pad", "Q", "", fit_bound_left, fit_bound_right);
  // _prf_graph_4pad->Fit("PRF_function_4pad", "Q", "", fit_bound_left, fit_bound_right);

  std::cout << "done" << std::endl;
  std::cout << "Process x histoes........................";

  for (auto i = 0; i < 4; ++i)
    ProfilePRF(_prf_histo_xscan[i], _prf_graph_xscan[i]);

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

  if (_prf_graph->GetFunction("PRF_function")) {
    _prf_function = _prf_graph->GetFunction("PRF_function");
    std::cout << "      PRF(x) = " << _prf_function->GetFormula()->GetExpFormula() << "  with ";
    for (auto i = 0; i < _prf_function->GetNpar(); ++i)
      std::cout << "  " << _prf_function->GetParameter(i) << ",";
    std::cout << "  Chi2/NDF " << _prf_graph->GetFunction("PRF_function")->GetChisquare()
              << "/" << _prf_graph->GetFunction("PRF_function")->GetNDF() << std::endl;
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

  if (_prev_iter_file)
    _prev_iter_file->Close();

  std::cout << "done" << std::endl;

  std::cout << "*************** Time consuming **************" << std::endl;
  std::cout << "Reading 3D array:\t"        << _sw_partial[5]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "Reconstruction:\t"        << _sw_partial[0]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "Analysis:      \t"        << _sw_partial[1]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "  Col loop:    \t"        << _sw_partial[2]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "  Fitters:     \t"        << _sw_partial[3]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "  Filling:     \t"        << _sw_partial[4]->CpuTime() * 1.e3 / (double)_eventList.size() << std::endl;
  return true;
}

//******************************************************************************
bool SpatialResolAna::ProfilePRF(const TH2F* PRF_h, TGraphErrors* gr) {
//******************************************************************************
  if (!PRF_h)
    return false;

  gr->Set(0);

  // Float_t threshold = 0.05 * PRF_h->GetMaximum();
  Float_t threshold = 0.;

  for (auto i = 1; i < PRF_h->GetXaxis()->GetNbins(); ++i) {

    TH1F* temp_h = (TH1F*)PRF_h->ProjectionY(Form("projections_bin_%i", i), i, i);
    if (!temp_h)
      return false;

    if (temp_h->GetMaximum() < threshold)
      continue;

    double x = PRF_h->GetXaxis()->GetBinCenter(i);
    double y = temp_h->GetBinCenter(temp_h->GetMaximumBin());

    gr->SetPoint(gr->GetN(), x, y);
    gr->SetPointError(gr->GetN()-1, 0, GetFWHM(temp_h)/2.);
  } // end of PRF histo profiling

  return true;
}

//******************************************************************************
bool SpatialResolAna::ProfilePRF_X(const TH2F* PRF_h, TGraphErrors* gr, TH1F* PRF_time_e) {
//******************************************************************************
  if (!PRF_h)
    return false;

  // Float_t threshold = 0.05 * PRF_h->GetMaximum();
  Float_t threshold = 1.;

  auto start = PRF_h->GetYaxis()->FindBin(0.);

  for (auto i = start+1; i < PRF_h->GetYaxis()->GetNbins(); ++i) {

    TH1F* temp_h = (TH1F*)PRF_h->ProjectionX(Form("projections_bin_%i", i), i, i);
    if (!temp_h)
      return false;

    if (temp_h->GetMaximum() < threshold)
      continue;

    for (auto j = temp_h->GetXaxis()->FindBin(0.); j < temp_h->GetXaxis()->GetNbins(); ++j)
      temp_h->SetBinContent(j, 0.);

    double y = PRF_h->GetYaxis()->GetBinCenter(i);
    double x = temp_h->GetBinCenter(temp_h->GetMaximumBin());

    gr->SetPoint(gr->GetN(), x, y);
    gr->SetPointError(gr->GetN()-1,  GetFWHM(temp_h)/2., 0.);
    PRF_time_e->Fill(y, GetFWHM(temp_h)/2);
  } // end of PRF histo profiling

  return true;
}

//******************************************************************************
Double_t SpatialResolAna::GetFWHM(const TH1F* h) {
//******************************************************************************
  auto mean = 0.;
  return GetFWHM(h, mean);
}

//******************************************************************************
Double_t SpatialResolAna::GetFWHM(const TH1F* h, Double_t& mean) {
//******************************************************************************
  if (!h->Integral())
    return -1.;

  mean = h->GetMean();
  auto max   = h->GetMaximum();
  auto start = h->GetBinLowEdge(h->FindFirstBinAbove(max/2));
  auto end   = h->GetBinLowEdge(h->FindLastBinAbove(max/2)) +
               h->GetBinWidth(h->FindLastBinAbove(max/2));

  return end - start;
}

//******************************************************************************
TF1* SpatialResolAna::InitializePRF(const TString name, bool shift) {
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

  TString formula;
  if (!shift)
    formula = "[0]*(1+[1]*x*x + [2] * x*x*x*x) / (1+[3]*x*x+[4]*x*x*x*x)";
  else
    formula = "[0]*(1+[1]*(x-[5])*(x-[5]) + [2] * (x-[5])*(x-[5])*(x-[5])*(x-[5])) / (1+[3]*(x-[5])*(x-[5])+[4]*(x-[5])*(x-[5])*(x-[5])*(x-[5]))";

  func  = new TF1(name,
    formula,
    prf_min, prf_max);
  func->SetParName(0, "Const");
  func->SetParName(1, "a2");
  func->SetParName(2, "a4");
  func->SetParName(3, "b2");
  func->SetParName(4, "b4");
  if (shift)
    func->SetParName(5, "shift");

  double co = 0.85436010;
  double a2 = -1919.1462;
  double a4 = 17575063.;
  double b2 = 627.88307;
  double b4 = 1.0339875e+09;
  double s = 0.;

  if (!shift) {
    func->SetParameters(co, a2, a4, b2, b4);
    return func;
  }

  func->SetParameters(co, a2, a4, b2, b4, s);
  return func;
}

bool SpatialResolAna::Draw() {
  std::cout << "Draw event " << _ev << std::endl;
  TCanvas c("c_SR", "c_SR", 800, 600);
  TCanvas c2("c_resid", "c_resid", 800, 0, 800, 600);
  TGraphErrors* gr = new TGraphErrors();
  TGraphErrors* gr_f = new TGraphErrors();
  TGraphErrors* gr_c = new TGraphErrors();
  TGraphErrors* gr_res_cl = new TGraphErrors();
  TGraphErrors* gr_res_av = new TGraphErrors();
  auto scale = 1.e3;
  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_clust_pos[colIt] == -999)
      continue;

    gr->SetPoint(gr->GetN(), scale*_x_av[colIt], scale*_clust_pos[colIt]);
    gr_f->SetPoint(gr_f->GetN(), scale*_x_av[colIt], scale*_track_pos[colIt]);

    gr_res_av->SetPoint(gr_res_av->GetN(),
                        scale*_x_av[colIt],
                        scale*(_clust_pos[colIt] - _track_pos[colIt]));
  }

  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_clust_pos[colIt] == -999)
      continue;

    gr_c->SetPoint(gr_c->GetN(), scale*_x[colIt], scale*_clust_pos[colIt]);

    double track_pos = 0.;
    if (_track_fit_func)
      track_pos = _track_fit_func->Eval(_x[colIt]);
    gr_res_cl->SetPoint(gr_res_cl->GetN(),
                        scale*_x[colIt],
                        scale*(_clust_pos[colIt] - track_pos));
  }

  c.cd();
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


  c2.cd();
  gPad->SetGrid();
  gr_res_cl->SetMarkerStyle(kPlus);
  gr_res_cl->GetYaxis()->SetTitle("Track Y, [mm]");
  gr_res_cl->GetXaxis()->SetTitle("Track X, [mm]");
  gr_res_cl->Draw("ap");

  gr_res_av->Draw("p same");
  c2.Update();


  c2.WaitPrimitive();
  return true;
}

//******************************************************************************
TCanvas* SpatialResolAna::DrawSelectionCan(const TRawEvent* event) {
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
  for (auto h:event->GetHits()){
    if(!h->GetQ()) continue;
    event3D->Fill(h->GetTime(),h->GetRow(),h->GetCol(),h->GetQ());
    MMsel->Fill(geom::y_pos[h->GetCol()],geom::x_pos[h->GetRow()],h->GetQ());
  }

  TCanvas *canv = new TCanvas("canv", "canv", 0., 0., 1400., 600.);
  canv->Divide(3,1);
  canv->cd(1);
  if (MM->Integral())
    MM->Draw("colz");
  //_prf_histo->Draw("COLZ");
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
}

//******************************************************************************
void SpatialResolAna::help(const std::string& name) {
//******************************************************************************
  AnalysisBase::help(name);
  std::cout << "   -t <iteration>       : iteration number" << std::endl;
  std::cout << "   --prev      <file>   : file from previous iteration" << std::endl;
  std::cout << std::endl;
}
