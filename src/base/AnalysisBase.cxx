#include <algorithm>
//#include <unistd.h>
#include <cstdlib>

#include "TROOT.h"

#include <GenericToolbox.h>

#include "AnalysisBase.hxx"

//******************************************************************************
AnalysisBase::AnalysisBase() :
  _clustering(nullptr),
  _file_in_name(""),
  _file_out_name(""),
  _param_file_name(""),
  _start_ID(-1),
  _end_ID(-1),
  _selected(0),
  _event(nullptr),
  _work_with_event_file(false),
  _file_in(nullptr),
  _file_out(nullptr),
  _chain(nullptr),
  _reconstruction(nullptr),
  _max_mult(6),
  _max_mean_mult(5),
  _cut_gap(true),
  _min_clusters(30),
  _verbose(1),
  _batch(false),
  _test_mode(false),
  _overwrite(false),
  _invert(false),
  _gaus_lorentz_PRF(false),
  _individual_column_PRF(false),
  _prf_free_centre(false),
  _do_linear_fit(false),
  _do_para_fit(false),
  _to_store_wf(true),
  _app(nullptr)
{
//******************************************************************************
  // CLI reader
  _clParser.setIsUnixGnuMode(true);
  _clParser.setIsFascist((true));

  _clParser.addOption("input_file", {"-i", "--input"}, "Input file name", 1);
  _clParser.addOption("output_file", {"-o", "--output"}, "Output file name", 1);

  _clParser.addOption("param_file", {"-p", "--param"}, "Parameter file name", 1);
  _clParser.addOption("verbosity", {"-v", "--verbose"}, "Verbosity level", 1);
  _clParser.addOption("start_id", {"--start"}, "Start event ID", 1);
  _clParser.addOption("end_id", {"--end"}, "End event ID", 1);

  _clParser.addTriggerOption("batch", {"-b"}, "Batch mode");
  _clParser.addTriggerOption("debug", {"-d", "--debug"}, "Debug mode");
  _clParser.addTriggerOption("overwrite", {"-r", "--overwrite"}, "Overwrite output file");

  _clParser.addTriggerOption("help", {"-h", "--help"}, "Print usage");
}

bool AnalysisBase::ReadCLI(int argc, char **argv) {
  _clParser.parseCmdLine(argc, argv);

  if (_clParser.isOptionTriggered("help")) {
    _clParser.getConfigSummary();
    exit(1);
  }

  setInputFile(_clParser.getOptionVal<TString>("input_file", "", 0));
  setOutputFile(_clParser.getOptionVal<TString>("output_file", "", 0));

  setParamFile(_clParser.getOptionVal<TString>("param_file", "", 0));
  setVerbosity(_clParser.getOptionVal<int>("verbosity", _verbose, 0));

  setStartID(_clParser.getOptionVal<int>("start_id", _start_ID, 0));
  setEndID(_clParser.getOptionVal<int>("end_id", _end_ID, 0));

  setBatchMode(_clParser.isOptionTriggered("batch"));
  setDebugMode(_clParser.isOptionTriggered("debug"));
  setOverwrite(_clParser.isOptionTriggered("overwrite"));

  if (!_batch)
    _app = new TApplication("app", &argc, argv);

  return true;
}

//******************************************************************************
bool AnalysisBase::Initialize() {
//******************************************************************************
  CL_col = new Clustering(0., 0);
  CL_diag = new Clustering(units::a45, 1);
  CL_2by1 = new Clustering(units::a2, 2);
  CL_3by1 = new Clustering(units::a3, 3);
  // CL_3by2 = new Clustering(units::a32, 2./3.);

  // Read parameter file
  if (!ReadParamFile()) {
    std::cerr << "ERROR! " << __func__ << "(). Parameter file is not read" << std::endl;
    exit(1);
  }

  if (_invert) {
    CL_2by1->angle = units::a2_inv;
    CL_3by1->angle = units::a3_inv;
    // CL_3by2->angle = units::a32_inv;
  }

  // read the first root file and decide
  // is it a raw file for the reconstruction or a file with TEvent
  TString tree_name = "";
  TFile* file;
  TString filename = _file_in_name;

  if (_file_in_name == "") {
    std::cerr << "ERROR. " << __func__ << "() No input file specified" << std::endl;
    exit(1);
  }
  if (_file_out_name == "") {
    std::cerr << "ERROR. " << __func__ << "() No output file specified" << std::endl;
    exit(1);
  }


  // extract the name of the input ROOT file
  // in case of list input take the first ROOT file
  if (!_file_in_name.Contains(".root")) {
    std::ifstream fList(_file_in_name.Data());
    if (!fList.good()) {
      std::cerr << "Can not read input " << _file_in_name << std::endl;
      exit(1);
    }
    std::string temp;
    getline(fList, temp);
    filename = static_cast<std::string>(temp);
  } // end of list input parse

  file = new TFile(filename.Data(), "READ");

  // find out which tree was send as an input
  // padAmpl[][][] or TRawEvent
  if ((TTree*)file->Get("tree")) {
    std::cout << "Raw data is using" << std::endl;
    tree_name = "tree";
    _work_with_event_file = false;
    ChainInputFiles(tree_name);
    if (!_chain) {
      std::cerr << "Error while chaining files" << std::endl;
      exit(1);
    }
    // Check the time binning. As we used both 510 and 511 time bins
    // the correct binning should be used for reading file
    TString branch_name = _chain->GetBranch("PadAmpl")->GetTitle();
    if (branch_name.Contains("[510]")) {
      _saclay_cosmics = true;
      _chain->SetBranchAddress("PadAmpl", _padAmpl_saclay);
    } else if (branch_name.Contains("[511]")) {
      _saclay_cosmics = false;
      _chain->SetBranchAddress("PadAmpl", _padAmpl);
    } else {
      std::cerr << "ERROR. AnalysisBase::Initialize()" << std::endl;
      std::cerr << "Time binning is inconsistent." << std::endl;
      std::cerr << "Read from file " << branch_name  << std::endl;
      exit(1);
    }
  } else if((TTree*)file->Get("event_tree")) {
    std::cout << "TRawEvent data is using" << std::endl;
    tree_name = "event_tree";
    ChainInputFiles(tree_name);
    _work_with_event_file = true;
    _chain->SetBranchAddress("Event", &_event);
  } else {
    std::cerr << "ERROR. AnalysisBase::Initialize. Unknown tree name" << std::endl;
    exit(1);
  }

  file->Close();

  std::cout << "Initializing analysis base...............";

  // setup the T2K style
  Int_t T2KstyleIndex = 2;
  // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
  TString localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  _t2kstyle = T2K().SetT2KStyle(T2KstyleIndex, localStyleName);

  gROOT->SetStyle(_t2kstyle->GetName());
  gROOT->ForceStyle();

  Long64_t N_events = _chain->GetEntries();
  for (auto i = 0; i < N_events; ++i)
    _eventList.push_back(i);

  // Open the output file
  if(_overwrite)
    _file_out = new TFile(_file_out_name.Data(), "RECREATE");
  else
    _file_out = new TFile(_file_out_name.Data(), "NEW");


  if (!_file_out->IsOpen()) {
    std::cerr << "ERROR. AnalysisBase::Initialize()" << std::endl;
    std::cerr << "File already exists or directory is not writable" << std::endl;
    std::cerr << "To prevent overwriting of the previous result the program will exit" << std::endl;
    exit(1);
  }


  if (_file_out)
    _file_out->cd();

  // Initialize histoes
  // * do it in your analysis *

  // Initialize selection
  // * do it in your analysis *

  std::cout << "done" << std::endl;

  return true;
}

//******************************************************************************
bool AnalysisBase::Loop() {
//******************************************************************************
  auto N_events = (int)_eventList.size();
  if (_test_mode)
    N_events = std::min((int)_eventList.size(), 100);

  if (_start_ID < 0)
    _start_ID = 0;
  if (_end_ID > 0) {
    _end_ID = std::min(_end_ID, N_events);
  } else
    _end_ID = N_events;

  _sw_event = new TStopwatch();

  _sw_partial[0] = new TStopwatch();
  _sw_partial[0]->Reset();
  _sw_partial[1] = new TStopwatch();
  _sw_partial[1]->Reset();
  _sw_partial[5] = new TStopwatch();
  _sw_partial[5]->Reset();

  if (_verbose >= v_progress) {
    std::cout << "Input file............................... " << _file_in_name << std::endl;
    std::cout << "Output file.............................. " << _file_out_name << std::endl;
    std::cout << "Processing" << std::endl;
    std::cout << "[                              ]   Nevents = " << _end_ID - _start_ID << "\r[";
    _sw_event->Start(false);
  }

  int denominator = 100;
  if (N_events < 100)
    denominator = N_events;

  // Event loop
  GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds(1);
  for (auto eventID = _start_ID; eventID < _end_ID; ++eventID) {
    if (_verbose >= v_event_number) {
      std::cout << "*************************************" << std::endl;
      std::cout << "Event " << eventID << std::endl;
      std::cout << "*************************************" << std::endl;
    }

    // Dump progress in command line
    if (_verbose == v_progress && (eventID%(N_events/denominator)) == 0)
      this->CL_progress_dump(eventID - _start_ID, _end_ID - _start_ID);

    _chain->GetEntry(_eventList[eventID]);

    _sw_partial[5]->Start(false);

    if (!_work_with_event_file) {
      // create TRawEvent from 3D array
      _event = std::make_shared<TRawEvent>(_eventList[eventID]);

      // Subtract the pedestal
      for (auto x = 0; x < geom::nPadx; ++x) {
        for (auto y = 0; y < geom::nPady; ++y) {
          auto hit = std::make_shared<THit>(x, y);
          // std::vector<int> adc;
          auto Qmax = -1;
          auto Tmax = -1;
          for (auto t = 0; t < geom::Nsamples; ++t) {
            // ommit last
            if (_saclay_cosmics && t == geom::Nsamples)
              continue;
            int q = _saclay_cosmics ?
                    _padAmpl_saclay[x][y][t] - 250 :
                    _padAmpl[x][y][t] - 250;

            hit->SetADC(t, q);
            if (q > Qmax) {
              Qmax = q;
              Tmax = t;
            }
            //
            /** REWEIGHT OF THE PAD*/
            // if (x == 5 && y == 16)
            //   _padAmpl[x][y][t] *= 0.95;
            /** */
          } // over time
          // hit->SetWF_v(adc);
          if (Qmax > 0){
            // compute FWHM
            int fwhm = 0;
            int width = 0;
            for (auto t = 0; t < geom::Nsamples; ++t) {
              if (hit->GetADC(t) > Qmax / 2)
                fwhm += 1;
              if (hit->GetADC(t) > 0)
                width += 1;
            }
            hit->SetFWHM(fwhm);
            hit->SetWidth(width);
            hit->SetQ(Qmax);
            hit->SetTime(Tmax);
            _event->AddHit(hit);
          } else {
            hit.reset();
          }
        } // over Y
      } // over X
    } // if 3D array input

    _sw_partial[5]->Stop();

    _store_event = false;

    _sw_partial[0]->Start(false);

    // copy event to a child class to be filled with reconstruction
    std::shared_ptr<TEvent> reco_event = std::make_shared<TEvent>(*_event);
    if (!_reconstruction->SelectEvent(reco_event)) {
      continue;
    }
    _sw_partial[0]->Stop();
    // do basic plotting
    auto c = std::make_unique<TCanvas>();
    if (!_batch) {
      c = DrawSelection(_event, reco_event);
      c->SetTitle(Form("Event %i", _event->GetID()));
      c->Draw();
      c->WaitPrimitive();
    }

    _sw_partial[1]->Start(false);
    ProcessEvent(reco_event);
    _sw_partial[1]->Stop();

    if (_store_event)
      ++_selected;
  } // end of event loop
//  std::cout << "time" << std::endl;
//   TODO move all time management to GT
//  std::cout << GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds(1) << std::endl;


  // if progress bar is active --> go to the next line
  if (_verbose == v_progress)
    std::cout << std::endl;

  return true;
}

//******************************************************************************
bool AnalysisBase::ProcessEvent(const std::shared_ptr<TEvent>& event) {
//******************************************************************************
  (void)event;
  throw std::logic_error("Event processing should be defined in your analysis");
}

//******************************************************************************
bool AnalysisBase::WriteOutput() {
//******************************************************************************
  //if(_test_mode) return true;
  if (!_file_out || !_file_out->IsOpen()){
    std::cout << "AnalysisBase::WriteOutput   _file_out is not Open!" << std::endl;
    return false;
  }

  std::cout << "Writing standard output..................";


  _file_out->cd();

  auto size = _output_vector.size();
  for (auto i = 0; i < size; ++i) {
    if (!_output_vector[i])
      std::cerr << "ERROR! " << __func__ << "()  output object pointer is nullptr" << std::endl;
    _output_vector[i]->Write();
  }

  _file_out->Close();

  std::cout << "done     " << "Write  " << size << " objects" << std::endl;

  return true;
}

//******************************************************************************
std::unique_ptr<TCanvas> AnalysisBase::DrawSelection(
    const std::shared_ptr<TRawEvent>& raw_event,
    const std::shared_ptr<TEvent>& reco_event
    ) {
//******************************************************************************
  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.05);
  TH2F    MM("MM","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TH2F    MMsel("MMsel","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TNtuple event3D("event3D", "event3D", "x:y:z:c");

//   all hits
  for (const auto& h : raw_event->GetHits()) {
    MM.Fill(h->GetCol(),h->GetRow(),h->GetQ());
  }

  // sel hits
  for (const auto& h : reco_event->GetUsedHits()){
    if(!h->GetQ()) continue;
    event3D.Fill((Float_t)h->GetTime(),(Float_t)h->GetRow(),(Float_t)h->GetCol(), (Float_t)h->GetQ());
    MMsel.Fill(h->GetCol(),h->GetRow(),h->GetQ());
  }

  for (auto x = 0; x < geom::nPadx; ++x) {
    for (auto y = 0; y < geom::nPady; ++y) {
      auto max = 0;
      for (auto t = 0; t < geom::Nsamples; ++t) {
        if (_padAmpl[x][y][t] > max) {
          max = _padAmpl[x][y][t];
        }
      } // over t
      if (max)
        MM.Fill(x, y, max);
    }
  }

  auto canv = std::make_unique<TCanvas>("canv", "canv", 0., 0., 1400., 600.);
  canv->Divide(3,1);
  canv->cd(1);
  MM.Draw("COLZ");
  canv->cd(2);
  MMsel.Draw("COLZ");

  canv->cd(3);
  event3D.Draw("x:y:z:c","","box2");
  auto htemp = (TH3F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetLimits(0,geom::nPadx);
  htemp->GetYaxis()->SetLimits(0,geom::nPady);
  htemp->GetZaxis()->SetLimits(0,500);
  htemp->SetTitle("");
  canv->Update();
  return canv;
}

//******************************************************************************
THitPtrVec AnalysisBase::GetRobustPadsInCluster(THitPtrVec col) {
//******************************************************************************
  std::vector<std::shared_ptr<THit>> result;
  // sort in charge decreasing order
  sort(col.begin(), col.end(), [](const std::shared_ptr<THit> & hit1,
                                  const std::shared_ptr<THit> & hit2) {
    return hit1->GetQ() > hit2->GetQ();
  });

  // leading pad
  auto col_id = col[0]->GetCol();
  auto row_id = col[0]->GetRow();
  // excluded from analysis the whole cluster if leading pad is near the broken pad
  for (const auto& broken : _broken_pads) {
    if (abs(col_id - broken.first) < 2 && abs(row_id - broken.second) < 2)
      return result;
  }

  for (const auto& pad : col) {
    auto q      = pad->GetQ();
    if (!q)
      continue;

    /** cross-talk candidate */
    // if (i > 0 &&
    //     pad->GetTime() - col[0]->GetTime() < 4 &&
    //     1.*q / col[0]->GetQ() < 0.08)
    //   continue;
    /** */

    // not more then 3 pads
    // if (i > 1)
    //   continue;

    // // WF with negative dt
    // if (pad->GetTime() - col[0]->GetTime() < -1)
    //   continue;

    // // avoid "suspicious" WF with small time difference in the 3rd pad
    // if (i > 1 && pad->GetTime() - col[0]->GetTime() < 5)
    //   continue;

    result.push_back(pad);

    // auto it_y   = pad->GetRow(_invert);
    // auto center_pad_y = geom::GetYpos(it_y, _invert);
  }

  return result;
}

//******************************************************************************
TClusterPtrVec AnalysisBase::GetRobustClusters(TClusterPtrVec & tr) {
//******************************************************************************
  TClusterPtrVec result;
  // sort clusters in increasing order
  sort(tr.begin(), tr.end(), [](TClusterPtr & cl1,
                                TClusterPtr & cl) {
                                        return  cl1->GetCharge() < cl->GetCharge();});

  // truncation cut
  /* NO TRUNCATION */
  auto frac = 1.00;
  auto i_max = (int)round(frac * (double)tr.size());
  result.reserve(i_max);
  for (auto i = 0; i < i_max; ++i) {
    result.push_back(std::move(tr[i]));
  }
  /* */

  /* truncate with prominence */
  // sort along the track
  // sort(tr.begin(), tr.end(), [](TCluster* cl1,
  //                               TCluster* cl){
  //                                 return  cl1->GetX() < cl->GetX();});
  // compute the prominence
  // auto prom_cut = 0.6;
  // for (uint i = 1; i < tr.size() - 1; ++i) {
  //   auto prom = 2.*tr[i]->GetCharge() / (tr[i-1]->GetCharge() + tr[i+1]->GetCharge());
  //   if (prom > prom_cut)
  //     result.push_back(tr[i]);
  // }
  /* */

  // BUG truncation with neighbours is not working with clusters
  // trancation + neibours
  // auto frac = 0.95;
  // std::vector<Int_t> bad_pads;
  // Int_t i_max = round(frac * tr.size());
  // for (uint i = i_max; i < tr.size(); ++i) {
  //   bad_pads.push_back(tr[i][0]->GetCol(_invert));
  // }
  // for (auto i = 0; i < i_max; ++i) {
  //   auto it_x = tr[i][0]->GetCol(_invert);
  //   if (find(bad_pads.begin(), bad_pads.end(), it_x+1) == bad_pads.end() ||
  //       find(bad_pads.begin(), bad_pads.end(), it_x-1) == bad_pads.end())
  //     result.push_back(tr[i]);
  // }

  // cut on the total charge in the cluster
  // auto q_cut = 2000;
  // for (auto col:tr) {
  //   auto total_q = accumulate(col.begin(), col.end(), 0,
  //                       [](const int& x, const THit* hit)
  //                       {return x + hit->GetQ();}
  //                       );
  //   if (total_q < q_cut)
  //     result.push_back(col);
  // }

  // sort by X for return
  sort(result.begin(), result.end(),
       [&](TClusterPtr & cl1, TClusterPtr & cl2){
          return cl1->GetX() < cl2->GetX();
        });
  return result;
}

//******************************************************************************
TClusterPtrVec AnalysisBase::ClusterTrack(const THitPtrVec &tr) const {
//******************************************************************************
  if (!_clustering) {
    std::cerr << "ERROR! AnalysisBase::ClusterTrack(). Clustering is not defined" << std::endl;
    exit(1);
  }
  TClusterPtrVec cluster_v;
  for (const auto& pad:tr) {
    auto col_id = pad->GetCol(_invert);
    auto row_id = pad->GetRow(_invert);

    // skip first and last row/column
    if (row_id == 0 || row_id == geom::GetNRow(_invert)-1 ||
        col_id == 0 || col_id == geom::GetNColumn(_invert)-1)
      continue;

    auto cons = _clustering->GetConstant(row_id, col_id);

    // search if the cluster is already considered
    TClusterPtrVec::iterator it;
    for (it = cluster_v.begin(); it < cluster_v.end(); ++it) {
      if (!(**it)[0]) {
        continue;
      }

      auto cluster_col = (**it)[0]->GetCol(_invert);
      auto cluster_row = (**it)[0]->GetRow(_invert);
      if (_clustering->GetConstant(cluster_row, cluster_col) == cons) {
        (*it)->AddHit(pad);
        (*it)->AddCharge(pad->GetQ());
        /** update X position */
        auto x_pad = geom::GetXposPad(pad, _invert, _clustering->angle);
        auto mult  = (*it)->GetSize();
        auto x_new = ((*it)->GetX() * ((Float_t)mult - 1) + x_pad) / (double)mult;
        (*it)->SetX((float_t)x_new);
        /** */

        break;
      }
    } // loop over track clusters
    // add new cluster
    if (it == cluster_v.end()) {
      auto first_cluster = std::make_unique<TCluster>(pad);
      first_cluster->SetX((float_t) geom::GetXposPad(pad, _invert, _clustering->angle));
      first_cluster->SetCharge(pad->GetQ());
      cluster_v.push_back(std::move(first_cluster));
    }
  } // over pads

  return cluster_v;
}

//******************************************************************************
//******************************************************************************
// ***** Functions below are utils: params reading, progress dump **************
//******************************************************************************
//******************************************************************************

//******************************************************************************
bool AnalysisBase::ReadParamFile() {
//******************************************************************************
  if (_param_file_name == ""){
    auto source = std::string(__FILE__);
    auto found = source.find_last_of('/');
    _param_file_name = source.substr(0, found) + "/../../params/default.ini";
  }
  std::cout << "*****************************************" << std::endl;
  std::cout << "Read parameters from " << _param_file_name << std::endl;
  std::ifstream cFile (_param_file_name);
  if (cFile.is_open()) {
    std::string line;
    while(getline(cFile, line)) {
      line.erase(std::remove_if(line.begin(),
                                line.end(),
                                isspace
                                ),
                line.end()
                );

      if(line[0] == '#' || line.empty())
        continue;
      auto delimiterPos = line.find('=');
      auto name = line.substr(0, delimiterPos);
      auto value = line.substr(delimiterPos + 1);
      // std::cout << name << " " << value << '\n';
      if (name == "cluster") {
        if (value == "column") {
          _clustering = CL_col;
          std::cout << "Column cluster is used" << std::endl;
        } else if (value == "diag") {
          _clustering = CL_diag;
          std::cout << "Diagonal cluster is used" << std::endl;
        } else if (value == "2by1") {
          _clustering = CL_2by1;
          std::cout << "2by1 cluster is used" << std::endl;
        } else if (value == "3by1") {
          _clustering = CL_3by1;
          std::cout << "3by1 cluster is used" << std::endl;
        // } else if (value == "3by2") {
        //   _clustering = CL_3by2;
        } else {
          std::cerr << "ERROR. Unknown clustering " << value << std::endl;
          return false;
        }
      } else if (name  == "invert") {
        if (value == "1") {
          _invert = true;
          std::cout << "Inverted geometry used" << std::endl;
        }
      } else if (name == "prf_shape") {
        if (value == "gaus_lorentz") {
          _gaus_lorentz_PRF = true;
          std::cout << "PRF is fit with Gaussian-Lorentzian" << std::endl;
        } else if (value == "pol4") {
          std::cout << "PRF is fit with 4th degree polinom" << std::endl;
        } else {
          std::cerr << "ERROR. Unknown PRF function " << value << std::endl;
          return false;
        }
      }  else if (name == "individual_prf") {
        if (value == "1") {
          _individual_column_PRF = true;
          std::cout << "Individual PRF for each column is used" << std::endl;
        }
      } else if (name == "prf_centre_freedom") {
        if (value == "1") {
          _prf_free_centre = true;
          std::cout << "PRF centre position is a free parameter of the fit" << std::endl;
        }
      } else if (name == "track_shape") {
        if (value == "parabola") {
          _do_para_fit = true;
          std::cout << "Parabola track fit is used" << std::endl;
        } else if (value == "linear") {
          _do_linear_fit = true;
          std::cout << "Linear track fit is used" << std::endl;
        } else if (value == "arc") {
          std::cout << "Arc track fit is used" << std::endl;
        } else {
          std::cerr << "ERROR. Unknown track shape " << value << std::endl;
          return false;
        }
      } else if (name == "max_mult") {
        _max_mult = TString(value).Atoi();
      } else if (name == "max_mean_mult") {
        _max_mean_mult = TString(value).Atof();
      } else if (name == "cut_gap") {
        if (value == "0")
          _cut_gap = false;
      } else if (name == "cluster_min") {
        _min_clusters = TString(value).Atoi();
      } else if (name == "max_phi") {
        _max_phi = (Float_t)TString(value).Atof();
      } else if (name == "max_theta") {
        _max_theta = (Float_t)TString(value).Atof();
      //switch to WF storage
      } else if (name == "to_store_wf") {
        if (value == "0") {
          _to_store_wf = false;
          std::cout << "WFs will NOT be stored" << std::endl;
        } else {
          std::cout << "WFs will be stored. Analysis will be slowed down" << std::endl;
        }
      } else if (name == "cross_talk") {
        if (value == "suppress") {
          _cross_talk_treat = suppress;
          std::cout << "Cross-talk will be suppressed" << std::endl;
        } else if (value == "cherry_pick") {
          _cross_talk_treat = cherry_pick;
          std::cout << "Cross-talk will be cherry-picked" << std::endl;
        } else {
          _cross_talk_treat = def;
          std::cout << "Cross-talk will not be treated" << std::endl;
        }
      } else if (name == "dead") {
        auto dead_pads = GenericToolbox::splitString(value, ";");
        if (!dead_pads.empty())
          std::cout << "Dead pads: ";
        for (const auto & pad : dead_pads) {
          auto coordinates = GenericToolbox::splitString(pad, ",");
          if (coordinates.size() != 2) {
            continue;
//            std::cerr << pad << std::endl;
//            throw std::logic_error("Wrong dead pad syntax");
          }
          _broken_pads.emplace_back(TString(coordinates[0]).Atoi(),
                                    TString(coordinates[1]).Atoi()
                                    );
          std::cout << _broken_pads.back().first << ", " << _broken_pads.back().second << "; ";
        }
        std::cout << std::endl;
      }
    }
  } else {
    return false;
  }
  std::cout << "*****************************************" << std::endl;
  return true;
}

//******************************************************************************
bool AnalysisBase::ChainInputFiles(const TString& tree_name) {
//******************************************************************************
  _chain = new TChain(tree_name);
  if (_file_in_name.Contains(".root")) {
    _chain->AddFile(_file_in_name);
  } else {
    std::ifstream fList(_file_in_name.Data());
    if (!fList.good()) {
      std::cerr << "Can not read input " << _file_in_name << std::endl;
      exit(1);
    }
    while (fList.good()) {
      std::string temp_filename;
      getline(fList, temp_filename);
      if (fList.eof()) break;
      _chain->AddFile(temp_filename.c_str());
    }
  }

  return true;
}

//******************************************************************************
void AnalysisBase::CL_progress_dump(int eventID, int N_events) {
//******************************************************************************
  auto mem = GenericToolbox::getProcessMemoryUsage();
  double CPUtime  = _sw_event->CpuTime();
  double REALtime = _sw_event->RealTime();
  int m, s;
  if (eventID) {
    int EET         = (int)((N_events - eventID) * REALtime / eventID);
    CPUtime *= 1.e3;  CPUtime /= eventID;
    m = EET / 60;
    s = EET % 60;
  }

  for (auto i = 0; i < 30; ++i)
    if (i < 30.*eventID/N_events) std::cout << "#";
    else std::cout << " ";
  std::cout << "]   Nevents = " << N_events << "\t" << round(1.*eventID/N_events * 100) << "%";
  std::cout << "\t Memory  " <<  mem / 1048576 << " " << "MB";
  std::cout << "\t Selected  " << _selected;
  if (eventID) {
    std::cout << "\t Av speed CPU " << CPUtime << " ms/event";
    std::cout << "\t EET real " << m << ":";
    if (s < 10)
      std::cout << "0";
    std::cout << s;
  }
  std::cout << "      \r[" << std::flush;
  _sw_event->Continue();
}
