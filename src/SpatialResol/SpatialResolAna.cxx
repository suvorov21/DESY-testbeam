/** @cond */
#include "TVector3.h"
#include "GenericToolbox.Root.h"
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
  _do_separate_pad_fit(false),
  _charge_uncertainty(true) {
//******************************************************************************
  _clParser.addOption("prev_file",
                      {"--prev"},
                      "File from previous iteration",
                      1
                      );

  _clParser.addOption("iter",
                      {"-t", "--iter"},
                      "Iteration number",
                      1
                      );

  _clParser.addTriggerOption("corr",
                             {"-c", "--corr"},
                             "Whether to apply correction"
                            );
}


//******************************************************************************
bool SpatialResolAna::ReadCLI(int argc, char **argv) {
//******************************************************************************
  AnalysisBase::ReadCLI(argc, argv);

  _prev_iter_name = _clParser.getOptionVal<TString>("prev_file", "", 0);
  _iteration = _clParser.getOptionVal<int>("iter", _iteration, 0);
  _correction = _clParser.isOptionTriggered("corr");

  return true;
}

//******************************************************************************
bool SpatialResolAna::Initialize() {
//******************************************************************************
  std::cout << "*****************************************" << std::endl;
  std::cout << "***   Spatial resolution analysis    ****" << std::endl;
  std::cout << "*****************************************" << std::endl;

  AnalysisBase::Initialize();

  std::cout << "Batch mode    :   " << _batch         << std::endl;
  std::cout << "Verbosity     :   " << _verbose       << std::endl;
  std::cout << "Iteration     :   " << _iteration     << std::endl;
  std::cout << "Debug         :   " << _test_mode     << std::endl;
  std::cout << "Correction    :   " << _correction     << std::endl;

  std::cout << "Initializing spatial resolution ana......";

  gErrorIgnoreLevel = kSysError;

  _uncertainty_vs_prf_gr_prev = nullptr;
  _uncertainty_vs_prf_histo   = nullptr;

  _prf_function = InitializePRF("PRF_function", _prf_free_centre);

  // Initialize graph for PRF profiling
  _prf_graph = new TGraphErrors();
  _prf_graph->SetName("PRF_graph");

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

    if (!_prf_function) {
      std::cerr << "ERROR. " << __func__ ;
      std::cerr << "  PRF function is not specified" << std::endl;
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
      for (auto rest = 0; rest < std::min(_clustering->n_pads, 3); ++rest) {
        _prf_function_arr[rest] = InitializePRF("PRF_function_tmp", true);

        TH2F* tmp = new TH2F("PRF_histo_tmp","", prf_bin, prf_min, prf_max, 150,0.,1.5);
        TString s = TString::Itoa(_clustering->n_pads, 10);
        TString r = TString::Itoa(rest, 10);
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
        TString s = TString::Itoa(_clustering->n_pads, 10);
        TString r = TString::Itoa(colId, 10);
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

    auto* oldtree = (TTree*)_prev_iter_file->Get("outtree");
    TH1F h("resol", "", resol_bin, resol_min, resol_max);
    oldtree->Project("resol", "residual");
    h.Fit("gaus", "Q", "");
    _uncertainty = (Float_t)h.GetFunction("gaus")->GetParameter(2);

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

  // Initialise histoes and graphs
  _prf_histo = new TH2F("PRF_histo","", prf_bin, prf_min, prf_max, 150,0.,1.5);

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

  _tree->Branch("max_mult",     &_m_max);
  _tree->Branch("mean_mult",    &_m_mean);

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

  _output_vector.push_back(_tree);

  // schedule the output for writing
  _output_vector.push_back(_prf_function);
  _output_vector.push_back(_prf_histo);
  _output_vector.push_back(_prf_graph);

  _uncertainty_vs_prf_gr = new TGraphErrors();
  _uncertainty_vs_prf_gr->SetName("uncertainty_vs_prf_gr");
  _output_vector.push_back(_uncertainty_vs_prf_gr);

  auto dir_resol = _file_out->mkdir("resol_column");
  _output_vector.push_back(dir_resol);


  auto dir_x_scan = _file_out->mkdir("x_scan");
  _output_vector.push_back(dir_x_scan);

  for (auto i = 0; i < prf_error_bins-1; ++i) {
    _uncertainty_prf_bins[i] = new TH1F(Form("error_prf_bin_%i", i),
      "", resol_bin, resol_min, resol_max);
    // _output_vector.push_back(_uncertainty_prf_bins[i]);
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
  _reconstruction = std::make_unique<DBSCANReconstruction>();
  _reconstruction->Initialize(_verbose);

  // Initialise track fitter
  TrackFitterBase::TrackShape shape = TrackFitterBase::arc;
  if (_do_linear_fit) {
    shape = TrackFitterBase::linear;
  } else if (_do_para_fit) {
    shape = TrackFitterBase::parabola;
  }

  _fitter = std::make_unique<TrackFitCern>(shape,
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
    _fitter->SetPRFarr(_prf_function_arr, _clustering->n_pads);
    _fitter->SetComplicatedPatternPRF(true);
  }

  if (_individual_column_PRF && _iteration) {
    _fitter->SetPRFarr(_prf_function_arr, 36);
    _fitter->SetIndividualPRF(true);
  }

  return true;
}

//******************************************************************************
bool SpatialResolAna::ProcessEvent(const std::shared_ptr<TEvent>& event) {
//******************************************************************************
  auto track_hits = event->GetUsedHits();
  if (track_hits.empty())
    return false;
  GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("sel");

  // reset tree values
  Reset((int)event->GetID());
// *****************************************************************************
// ************** See steps documentation in the README.md file ****************
// *******************  STEP 1 *************************************************

  TClusterPtrVec clusters;
  clusters = ClusterTrack(track_hits);

  if (_verbose >= v_analysis_steps)
    std::cout << "Clusterization done " << clusters.size() <<  std::endl;

  // selection
  bool sel = sel::CrossingTrackSelection(clusters,
                                   _max_mult,
                                   _max_mean_mult,
                                   _cut_gap,
                                   _max_phi,
                                   _max_theta,
                                   _broken_pads,
                                   _invert,
                                   _verbose
                                   );
  _sel_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("sel");
  if (!sel)
    return false;

  // TODO prevent multiple fitter call
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
         [&](const TClusterPtr & cl1, const TClusterPtr & cl2){
            return cl1->GetX() < cl2->GetX();
          });
    clusters.erase(clusters.begin(), clusters.begin() + 1);
    clusters.erase(clusters.end() - 1);
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
  GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("column");
// *******************  STEP 2 *************************************************

  if (_verbose >= v_analysis_steps)
    std::cout << "start cluster fit" << std::endl;

  std::vector<int> QsegmentS;

  _m_mean = (Float_t)GenericToolbox::getAverage(robust_clusters, [](auto cluster) {return cluster->GetSize();});
  _m_max = (int)(*(std::max_element(robust_clusters.begin(),
                                    robust_clusters.end(),
                                    [](const auto& a, const auto& b){
                                      return a->GetSize() > b->GetSize();
                                    }
                                    )))->GetSize();

  for (uint clusterId = 0; clusterId < robust_clusters.size(); ++clusterId) {
    if (!(*robust_clusters[clusterId])[0])
      continue;

    // loop over rows
    auto robust_pads = GetRobustPadsInCluster(robust_clusters[clusterId]->GetHits());
    _multiplicity[clusterId] = (int)robust_pads.size();

    auto pad_id = -1;
    for (auto & pad:robust_pads) {
      ++pad_id;
      if (!pad)
        continue;

      // treat cross-talk
      // not for the first pad
      if (pad != robust_pads[0] && _cross_talk_treat != def)
        TreatCrossTalk(pad, robust_pads, pad_id);

      FillPadOutput(pad, (int)clusterId, pad_id);
    } // loop over pads

    _clust_pos[clusterId] /= (Float_t)_charge[clusterId];
    _x[clusterId] = robust_clusters[clusterId]->GetX();

    double CoC =  _clust_pos[clusterId];

    if (_multiplicity[clusterId] > 1 && _iteration > 0) {
      _clust_pos[clusterId] = (Float_t)_fitter->FitCluster(robust_pads,
                                                 _charge[clusterId],
                                                 _clust_pos[clusterId]
                                                 );
    }

    robust_clusters[clusterId]->SetY(_clust_pos[clusterId]);
    robust_clusters[clusterId]->SetCharge(_charge[clusterId]);

    if (_verbose >= v_fit_details) {
      for (const auto& pad:*robust_clusters[clusterId]) {
        std::cout << pad->GetRow(_invert) <<  " : " << pad->GetCol(_invert) << "\t";
      }
      std::cout << std::endl;
    }

    if (_verbose >= v_fit_details)
      std::cout << "X:CoC:Fit\t" << robust_clusters[clusterId]->GetX() << "\t" << CoC << "\t" << robust_clusters[clusterId]->GetY() << std::endl;

    if (robust_clusters[clusterId]->GetSize() == 1) {
      robust_clusters[clusterId]->SetYE(0.002);
    } else {
      if (_iteration)
        robust_clusters[clusterId]->SetYE(_uncertainty);
      else
        robust_clusters[clusterId]->SetYE(0.001);
    }

    if (_verbose >= v_residuals)
      std::cout << "Cluster pos " << _clust_pos[clusterId] << "\t";

    if (robust_clusters[clusterId]->GetCharge()) QsegmentS.push_back(robust_clusters[clusterId]->GetCharge());

  } // loop over clusters

  if (_verbose >= v_analysis_steps) {
    std::cout << "Loop over columns done" << std::endl;
  }

  _dEdx = ComputedEdx(QsegmentS);

  _column_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("column");
  GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("fitter");

//******************** STEP 3 **************************************************
  auto fit = _fitter->FitTrack(robust_clusters);
  _track_fit_func = fit;

  if (!fit)
    return false;

  _quality = (Float_t)fit->GetChisquare() / (Float_t)fit->GetNDF();

  if (_verbose >= v_analysis_steps)
    std::cout << "Track fit done" << std::endl;

  TString func = fit->GetName();

  // in case of arc fitting fill the momentum histo
  if (!_do_linear_fit && !_do_para_fit) {
    _mom = Float_t(1./fit->GetParameter(0) * units::B * units::clight / 1.e9);
    _sin_alpha  = (Float_t)fit->GetParameter(1);
    _offset     = (Float_t)fit->GetParameter(2);
  }

  if (_do_para_fit) {
    // pol2 := [p0]+[p1]*x+[p2]*pow(x,2)

    double x1 = robust_clusters[0]->GetX();
    TVector3 start(x1, fit->Eval(x1), 0.);
    double x2 = (*(robust_clusters.end()-1))->GetX();
    TVector3 end(x2, fit->Eval(x2), 0.);
    TLine_att line(start, end);

    double a = fit->GetParameter(2);
    double b = fit->GetParameter(1);

    double x0 = - b * (x2-x1) + a*x2*x2 + b*x2 - a*x1*x1 - b*x1 ;
    x0 /= 2*a*(x2-x1);
    TVector3 max_sag(x0, fit->Eval(x0), 0.);
    double sag = line.GetDistVec(max_sag).Mag();
    double L = (end - start).Mag();

    if (abs(sag) > 1e-10) {
      _mom = (Float_t)(L*L /8/sag + sag / 2);
      _mom *= (Float_t)(units::B * units::clight / 1.e9);
      if (a < 0)
        _mom *= -1;
    }

    _sin_alpha = (Float_t)TMath::Sin(TMath::ATan(fit->Derivative(x1)));
    _offset    = (Float_t)start.Y();

    if (_verbose > v_analysis_steps) {
      std::cout << "start:end:max\t" << start.X() << ", " << start.Y();
      std::cout << "\t" << end.X() << ", " << end.Y();
      std::cout << "\t" << max_sag.X() << ", " << max_sag.Y() << std::endl;
      std::cout << "Line at x0\t" << line.EvalX((Float_t)max_sag.X()).Y() << std::endl;
      std::cout << "Sagitta:\t" << sag << "\tLength:\t" << L << std::endl;
      std::cout << "mom:\t" << _mom << std::endl;
    }
  }

//****************** STEP 4 ****************************************************

  std::array<std::shared_ptr<TF1>, Nclusters> fit1{nullptr};
  for (int i = 0; i < Nclusters; ++i) {
    if (!_correction)
      fit1[i] = fit;
    else {
      fit1[i] = _fitter->FitTrack(robust_clusters, i);
    }
  }

  _fitters_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("fitter");
  GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("filling");

  // second loop over columns
  for (uint clusterId = 0; clusterId < robust_clusters.size(); ++clusterId) {
    if (!robust_clusters[clusterId]) continue;

//****************** STEP 5 ****************************************************
    FillSR(robust_clusters[clusterId], clusterId, fit, fit1);
// ************ STEP 6 *********************************************************
    // don't fill PRF for multiplicity of 1
    if (robust_clusters[clusterId]->GetSize() == 1)
      continue;
    // Fill PRF
    // TODO cache the result of GetRobustPadsInCluster?
    auto robust_pads = GetRobustPadsInCluster(robust_clusters[clusterId]->GetHits());
    int padId = 0;
    for (const auto& pad:robust_pads)
      FillPRF(pad, padId, clusterId, fit);
  } // loop over columns

  _filling_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("filling");

  if (!_batch)
    Draw();

  _tree->Fill();

  if (_store_event)
    _passed_events.push_back((Int_t)event->GetID());

  return true;
}

//******************************************************************************
void SpatialResolAna::TreatCrossTalk(const THitPtr& pad,
                                     THitPtrVec& robust_pads,
                                     int& pad_id) {
  //******************************************************************************
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
  } // cross-talk cherry-picking
}

//******************************************************************************
void SpatialResolAna::FillPadOutput(const THitPtr &pad,
                                    const int& clusterId,
                                    const int& padId) {
//******************************************************************************
  _clust_pos[clusterId] += (Float_t)pad->GetQ() * (Float_t)geom::GetYposPad(pad,
                                                                            _invert,
                                                                            _clustering->angle);
  _charge[clusterId] += pad->GetQ();
  _pad_charge[clusterId][padId] = pad->GetQ();
  _pad_time[clusterId][padId] = pad->GetTime();
  _pad_x[clusterId][padId] = pad->GetCol(_invert);
  _pad_y[clusterId][padId] = pad->GetRow(_invert);
  _wf_width[clusterId][padId] = pad->GetWidth();
  _wf_fwhm[clusterId][padId] = pad->GetFWHM();

  if (_to_store_wf) {
    for (int tz = 0; tz < geom::Nsamples; ++tz) {
      _pad_wf_q[clusterId][padId][tz] = pad->GetADC(tz);
    }
  }
}

//******************************************************************************
Float_t SpatialResolAna::ComputedEdx(std::vector<int>& QsegmentS) {
//******************************************************************************
  Float_t alpha = 0.7;
  sort(QsegmentS.begin(), QsegmentS.end());
  double totQ = 0.;
  auto i_max = (int)round(alpha * (double)QsegmentS.size());
  for (auto i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i)
    totQ += QsegmentS[i];

  return Float_t(totQ / (alpha * (Float_t)QsegmentS.size()));

}

//******************************************************************************
void SpatialResolAna::FillSR(const TClusterPtr& cluster,
                             const uint& clusterId,
                             const std::shared_ptr<TF1>& fit,
                             const std::array<std::shared_ptr<TF1>, Nclusters>& fit1
                             ) {
//******************************************************************************
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
  _track_pos[clusterId] = (Float_t)track_fit_y;
  _residual[clusterId] = (Float_t)(_clust_pos[clusterId] - track_fit_y);
  _residual_corr[clusterId] = (Float_t)(_clust_pos[clusterId] - track_fit_y1);
}

//******************************************************************************
void SpatialResolAna::FillPRF(const THitPtr& pad,
                              int& padId,
                              const uint& clusterId,
                              const std::shared_ptr<TF1>& fit) {
//******************************************************************************
  /// WARNING only 10 first pads are used
  if (!pad || padId > 9)
    return;

  auto q    = pad->GetQ();
  auto time = pad->GetTime();

  double x = geom::GetXposPad(pad, _invert, _clustering->angle);
  double center_pad_y = geom::GetYposPad(pad,
                                         _invert,
                                         _clustering->angle
  );

  double track_fit_y_pad    = fit->Eval(x);

  // Pad characteristics
  // time, position wrt track, charge, col/row
  _time[clusterId][padId]  = time;
  _dx[clusterId][padId]    = (Float_t)(center_pad_y - track_fit_y_pad);
  _qfrac[clusterId][padId] = (Float_t)q / (Float_t)_charge[clusterId];

  if (_verbose >= v_prf) {
    std::cout << "PRF fill\t" << _dx[clusterId][padId];
    std::cout << "\t" << _qfrac[clusterId][padId];
    std::cout << "\t" << _time[clusterId][padId] - _time[clusterId][0] << std::endl;
  }

  // fill PRF
  _prf_histo->Fill( _dx[clusterId][padId],
                   _qfrac[clusterId][padId]
  );

  if (padId > 0)
    _prf_time->Fill(_dx[clusterId][padId],
                    _time[clusterId][padId] - _time[clusterId][0]);

  ++padId;
}

//******************************************************************************
void SpatialResolAna::Reset(int id) {
//******************************************************************************
  _ev = id;
  _quality    = -999.;
  _mom        = -999.;
  _sin_alpha  = -999.;
  _offset     = -999.;
  _rob_clusters = -999;

  for (auto colId = 0; colId < Nclusters; ++colId) {
    _multiplicity[colId]  = -999;
    _charge[colId]        = 0;
    _residual[colId]      = -999;
    _residual_corr[colId] = -999;
    _clust_pos[colId]     = 0;
    _track_pos[colId]     = -999;
    _x[colId]             = -999;
    _x_av[colId]          = -999;
    _dEdx               = -999;

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
        _pad_wf_q[colId][padId][nSmp] = -999;
      }
    }
  }
}

//******************************************************************************
bool SpatialResolAna::WriteOutput() {
//******************************************************************************
  if (!_file_out)
    return true;

  std::cout << "PRF profiling............................";

  ProfilePRF(_prf_histo, _prf_graph);

  std::cout << "done" << std::endl;
  std::cout << "PRF fit..................................";

  // MAGIC
  for (auto i = 0; i < 3; ++i)
    _prf_graph->Fit("PRF_function", "Q", "", fit_bound_left, fit_bound_right);


  std::cout << "done" << std::endl;

  std::cout << "done" << std::endl;
  std::cout << "Process histoes..........................";

  if (_do_separate_pad_fit && _iteration) {
    for (auto res : _uncertainty_prf_bins) {
      Double_t mean = res->GetMean();
      Double_t sigma = 0.5*GenericToolbox::getFWHM(res);

      res->Fit("gaus", "Q", "", mean - 4*sigma, mean + 4*sigma);

      TF1* func = res->GetFunction("gaus");

      if (!func) {
        std::cout << "WARNING. SpatialResolAna::WriteOutput(). Residual fit fail" << std::endl;
        exit(1);
      }

      Double_t sigma_e;

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
    TH1F sr_h("h", "", resol_bin, resol_min, resol_max);
    _tree->Project("h", "residual");
    sr_h.Fit("gaus", "Q");
    std::cout << "Spatial resolution\t" << sr_h.GetFunction("gaus")->GetParameter(2) * 1e6 << " um" << std::endl;
    TH1F e_h("h_e", "", 1000, 0., 10000);
    _tree->Project("h_e", "dEdx");
    e_h.Fit("gaus", "Q");
    auto resol = e_h.GetFunction("gaus")->GetParameter(2) / e_h.GetFunction("gaus")->GetParameter(1);
    std::cout << "dE/dx             \t" << resol * 100 << "%" << std::endl;
  }

  // Write objects
  AnalysisBase::WriteOutput();

  std::cout << "Writing spatial analisis output..........";

  auto file = new TFile(_file_out_name.Data(), "UPDATE");
  // write
  auto tree = new TTree("EventTree", "");
  Int_t var = 0;
  tree->Branch("PassedEvents",     &var);
  for (int _passed_event : _passed_events) {
    var = _passed_event;
    tree->Fill();
  }
  tree->Write("", TObject::kOverwrite);
  file->Close();

  if (_prev_iter_file)
    _prev_iter_file->Close();

  std::cout << "done" << std::endl;

  std::cout << "************** Time consumption *************" << std::endl;
  std::cout << "Reading 3D array:\t" << (double)_read_time / 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "Reconstruction:  \t" << (double)_reco_time / 1.e3 / (double)_eventList.size() << std::endl;
  std::cout << "Analysis:        \t" << (double)_ana_time / 1.e3 / (double)_reconstructed << std::endl;
  std::cout << "  Selection:     \t" << (double)_sel_time / 1.e3 / (double)_reconstructed << std::endl;
  std::cout << "  Col loop:      \t" << (double)_column_time / 1.e3 / (double)_selected << std::endl;
  std::cout << "  Fitters:       \t" << (double)_fitters_time / 1.e3 / (double)_selected << std::endl;
  std::cout << "  Filling:       \t" << (double)_filling_time / 1.e3 / (double)_selected << std::endl;
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
    gr->SetPointError(gr->GetN()-1, 0, GenericToolbox::getFWHM(temp_h)/2.);
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
    gr->SetPointError(gr->GetN()-1,  GenericToolbox::getFWHM(temp_h)/2., 0.);
    PRF_time_e->Fill(y, GenericToolbox::getFWHM(temp_h)/2);
  } // end of PRF histo profiling

  return true;
}

//******************************************************************************
TF1* SpatialResolAna::InitializePRF(const TString& name, bool shift) {
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

//******************************************************************************
bool SpatialResolAna::Draw() {
//******************************************************************************
  std::cout << "Draw event " << _ev << std::endl;
  TCanvas c("c_SR", "c_SR", 0, 600, 800, 600);
  TCanvas c2("c_resid", "c_resid", 800, 600, 800, 600);
  auto* gr = new TGraphErrors();
  auto* gr_f = new TGraphErrors();
  auto* gr_c = new TGraphErrors();
  auto* gr_res_cl = new TGraphErrors();
  auto* gr_res_av = new TGraphErrors();
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
  gr_c->SetTitle("Event " + TString::Itoa(_ev, 10));
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
