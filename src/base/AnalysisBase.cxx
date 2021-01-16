#include <algorithm>
#include <unistd.h>
#define GetCurrentDir getcwd

#include "TROOT.h"

#include "AnalysisBase.hxx"

//******************************************************************************
AnalysisBase::AnalysisBase(int argc, char** argv) :
  _file_in_name(""),
  _file_out_name(""),
  _param_file_name(""),
  _start_ID(-1),
  _end_ID(-1),
  _selected(0),
  _event(NULL),
  _store_event_tree(false),
  _work_with_event_file(false),
  _file_in(NULL),
  _file_out(NULL),
  _chain(NULL),
  _Prev_iter_name(TString("")),
  _iteration(0),
  _reconstruction(NULL),
  _max_mult(6),
  _cut_gap(true),
  _min_clusters(30),
  _verbose(1),
  _batch(false),
  _test_mode(false),
  _overwrite(false),
  _invert(false),
  _gaus_lorentz_PRF(false),
  _do_linear_fit(false),
  _do_para_fit(false),
  _app(NULL),
  _useCern(false),
  _to_store_wf(true)
{
//******************************************************************************

  // TODO redefine CLI.
  // now written in an ugly way
  // prevent copy-past between the daughter-parent classes
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

    {"param",           required_argument, 0,   0},         // 8

    {"prev",            required_argument, 0,   0},         // 9

    {"help",            no_argument,    0,    'h'},         // 10

    {0,                 0,              0,      0}
  };

  int index;

  // read CLI
  for (;;) {
    int c = getopt_long(argc, argv, "i:o:bv:drhst:cp:", longopts, &index);
    if (c < 0) break;
    switch (c) {
      case 0  :
        if (index == 6) _start_ID         =  atoi(optarg);
        if (index == 7) _end_ID           =  atoi(optarg);
        if (index == 8) _param_file_name  = optarg;
        if (index == 9) _Prev_iter_name   = optarg;
        if (index == 10) help(argv[0]);
        break;
      case 'i' : _file_in_name     = optarg;       break;
      case 'o' : _file_out_name    = optarg;       break;
      case 't' : _iteration        = atoi(optarg); break;
      case 'b' : _batch            = true;         break;
      case 'v' : _verbose          = atoi(optarg); break;
      case 'd' : _test_mode        = true;         break;
      case 'p' : _param_file_name = optarg;        break;
      case 'r' :
        _overwrite        = true;
        std::cout << "Output will be overwritten" << std::endl;
        break;
      case 's' :
        _store_event_tree = true;
        std::cout << "Tree with TEvents will be written" << std::endl;
        break;
      case 'h' : help(argv[0]);                    break;
      //case '?' : help(argv[0]);
    }
  }

  if (_file_in_name == "") {
    std::cerr << "ERROR. AnalysisBase::AnalysisBase. No input file specified" << std::endl;
    exit(1);
  }

  if (_iteration == -1) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna().";
    std::cout << " Iteration should be defined as a input param" << std::endl;
    exit(1);
  }

  if (!_batch)
    _app = new TApplication("app", &argc, argv);
}

//******************************************************************************
bool AnalysisBase::Initialize() {
//******************************************************************************

  // WARNING
  // A very dirty adoptation of angles
  CL_col = new Clustering(0., 0);
  CL_diag = new Clustering(units::a45, 1);
  CL_2by1 = new Clustering(units::a2, 2);
  CL_3by1 = new Clustering(units::a3, 3);


  // Read parameter file
  if (!ReadParamFile()) {
    std::cerr << "ERROR! AnalysisBase::Initialize(). Parameter file is not read" << std::endl;
    exit(1);
  }

  if (_invert) {
    CL_2by1->angle = units::a2_inv;
    CL_3by1->angle = units::a3_inv;
  }

  // read the first root file and decide
  // is it a raw file for the reconstruction or a file with TEvent
  TString tree_name = "";
  TFile* file;
  TString filename = _file_in_name;

  // in case of list input
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
  // padAmpl or TEvent

  if ((TTree*)file->Get("tree")) {
    std::cout << "Raw data is using" << std::endl;
    tree_name = "tree";
    _work_with_event_file = false;
  } else if((TTree*)file->Get("event_tree")) {
    std::cout << "TEvent data is using" << std::endl;
    tree_name = "event_tree";
    _work_with_event_file = true;
  } else if ((TTree*)file->Get("padData")) {
    _useCern = true;
    // tree_name = "padData";
    // std::cout << "Running over CERN data" << std::endl;
    // _tgeom= (TTree*)file->Get("femGeomTree");
    // //std::vector<int> *iPad(0);
    // //std::vector<int> *jPad(0);
    // _tgeom->SetBranchAddress("jPad", &_jPad );
    // _tgeom->SetBranchAddress("iPad", &_iPad );
    // _tgeom->GetEntry(0); // put into memory geometry info
  } else {
    std::cerr << "ERROR. AnalysisBase::Initialize. Unknown tree name" << std::endl;
    exit(1);
  }

  file->Close();

  if (_work_with_event_file && _store_event_tree) {
    std::cerr << "ERROR. AnalysisBase::Initialize. Prohibited to generate TEvent over TEvent. Exit" << std::endl;
    exit(1);
  }

  std::cout << "Initializing analysis base...............";
  // read and chain input files
  _chain = new TChain(tree_name);
  TString first_file_name = "";

  if (_file_in_name.Contains(".root")) {
    _chain->AddFile(_file_in_name);
    first_file_name = _file_in_name;
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
      if (first_file_name.CompareTo("") == 0)
        first_file_name = temp_filename;
    }
  }

  if (_useCern) {
    std::cerr << "CERN data format is deprecated" << std::endl;
    exit(1);
    // _chain->SetBranchAddress("PadphysChannels", &_listOfChannels );
    // _chain->SetBranchAddress("PadADCvsTime"   , &_listOfSamples );

  } else if (!_work_with_event_file) {
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
  } else
    _chain->SetBranchAddress("Event", &_event);

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

  Int_t N_events = _chain->GetEntries();
  for (auto i = 0; i < N_events; ++i)
    _EventList.push_back(i);

  // Open the output file
  //if (!_test_mode){
  if (1){

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
  }

  // in case we want to store TEvent in the file
  if (_store_event_tree) {
    // take a name from the input and dir from output
    Ssiz_t slash_pos = 0;
    while (_file_out_name.Index("/", 1, slash_pos+1, TString::kExact) != -1)
      slash_pos = _file_out_name.Index("/", 1, slash_pos+1, TString::kExact);
    TString file_dir  = _file_out_name(0, slash_pos+1);
    slash_pos = 0;
    while (first_file_name.Index("/", 1, slash_pos+1, TString::kExact) != -1)
      slash_pos = first_file_name.Index("/", 1, slash_pos+1, TString::kExact);
    TString file_name = first_file_name(slash_pos+1, first_file_name.Length());

    _event_file = new TFile((file_dir + file_name).Data(), "RECREATE");
    _event_tree = new TTree("event_tree", "");
    _event_tree->Branch("Event",    &_event,  32000,  0);
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
bool AnalysisBase::Loop(std::vector<Int_t> EventList) {
//******************************************************************************
  auto N_events = static_cast<Int_t>(EventList.size());
  if (_test_mode)
    N_events = std::min(static_cast<Int_t>(EventList.size()), 100);

  if (_start_ID < 0)
    _start_ID = 0;
  if (_end_ID > 0)
    N_events = _end_ID;

  _sw_event = new TStopwatch();

  _sw_partial[0] = new TStopwatch();
  _sw_partial[0]->Reset();
  _sw_partial[1] = new TStopwatch();
  _sw_partial[1]->Reset();

  if (_verbose >= v_progress) {
    std::cout << "Input file..............................." << _file_in_name << std::endl;
    std::cout << "output file.............................." << _file_out_name << std::endl;
    std::cout << "Processing" << std::endl;
    std::cout << "[                              ]   Nevents = " << N_events << "\r[";
    _sw_event->Start(0);
  }

  int denimonator = 100;
  if (N_events < 100)
    denimonator = N_events;
  for (auto eventID = _start_ID; eventID < N_events; ++eventID) {
    if (_verbose >= v_event_number) {
      std::cout << "*************************************" << std::endl;
      std::cout << "Event " << eventID << std::endl;
      std::cout << "*************************************" << std::endl;
    }

    if (_verbose == v_progress && (eventID%(N_events/denimonator)) == 0)
      this->CL_progress_dump(eventID - _start_ID, N_events - _start_ID);

    _chain->GetEntry(EventList[eventID]);

    if (_useCern) {
      std::cerr << "CERN data format is deprecated" << std::endl;
      exit(1);
      // memset(_padAmpl, 0, geom::nPadx * geom::nPady * geom::Nsamples * (sizeof(Int_t)));
      // for (uint ic=0; ic< _listOfChannels->size(); ic++){
      //   int chan= (*_listOfChannels)[ic];
      //   if ((*_jPad)[chan] >= geom::nPadx ||
      //         (*_iPad)[chan] >= geom::nPady)
      //     continue;
      //   for (uint it = 0; it < (*_listOfSamples)[ic].size(); it++){
      //     if (it >= geom::Nsamples)
      //       continue;
      //     int adc = (*_listOfSamples)[ic][it];
      //     _padAmpl[(*_jPad)[chan]][(*_iPad)[chan]][it] = adc;
      //   }
      // }
    } else {
      // Subtract the pedestal
      for (auto x = 0; x < geom::nPadx; ++x) {
        for (auto y = 0; y < geom::nPady; ++y) {
          for (auto t = 0; t < geom::Nsamples; ++t) {
            // ommit last
            if (_saclay_cosmics && t == geom::Nsamples)
              continue;
            int q = 0;
            if (_saclay_cosmics) {
              q = _padAmpl_saclay[x][y][t] - 250;
            } else {
              q = _padAmpl[x][y][t] - 250;
            }

            //_padAmpl[x][y][t] = q < 0 ? 0 : q;
            _padAmpl[x][y][t] = q < -249 ? 0 : q;

            /** REWEIGHT OF THE PAD*/
            // if (x == 5 && y == 16)
            //   _padAmpl[x][y][t] *= 0.95;
            /** */

          }
        }
      }
    }

    _store_event = false;

    _sw_partial[0]->Start(false);

    if (!_work_with_event_file) {
      if (_event && !_store_event_tree)
        delete _event;
      _event = new TEvent(EventList[eventID]);

      if (!_reconstruction->SelectEvent(_padAmpl, _event))
        continue;
    }
    // else _event->SetID(EventList[eventID]);

    _sw_partial[0]->Stop();
    _sw_partial[1]->Start(false);
    ProcessEvent(_event);
    _sw_partial[1]->Stop();

    if (_store_event) {
      ++_selected;
      if (_store_event_tree)
        _event_tree->Fill();
    }

    if (!_store_event_tree && !_work_with_event_file) {
      delete _event;
      _event = NULL;
    }
  } // end of event loop

  if (_verbose == v_progress)
    std::cout << std::endl;

  return true;
}

//******************************************************************************
bool AnalysisBase::ProcessEvent(const TEvent* event) {
//******************************************************************************
  (void)event;
  std::cerr << "EROOR. AnalysisBase::ProcessEvent(). Event processing should be defined in your analysis" << std::endl;
  exit(1);
  return true;
}

//******************************************************************************
bool AnalysisBase::WriteOutput() {
//******************************************************************************
  //if(_test_mode) return true;
  if (!_file_out->IsOpen()){
    std::cout << "AnalysisBase::WriteOutput   _file_out is not Open!" << std::endl;
    return false;
  }

  // Write the TEvents in the file
  if (_store_event_tree) {
    _event_file->cd();
    _event_tree->Write("", TObject::kOverwrite);
    std::cout << "Wrote TEvent events into " << _event_file->GetName() << std::endl;
    _event_file->Close();
  }

  std::cout << "Writing standard output..................";


  _file_out->cd();

  auto size = static_cast<int>(_output_vector.size());
  for (auto i = 0; i < size; ++i) {
    if (!_output_vector[i])
      std::cerr << "ERROR! AnalysisBase::WriteOutput()  output object pointer is NULL" << std::endl;
    _output_vector[i]->Write();
  }

  _file_out->Close();

  std::cout << "done     " << "Write  " << size << " objects" << std::endl;

  return true;
}

// TODO
// make the inheritance possible
// e.g. draw events here but also draw some analysi specific stuff in the analysis
//******************************************************************************
void AnalysisBase::DrawSelection(const TEvent *event, int trkID){
//******************************************************************************
  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.05);
  TH2F    *MM      = new TH2F("MM","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TH2F    *MMsel   = new TH2F("MMsel","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");

  // all hits
  //for(auto h:event->GetHits()){
  //  MM->Fill(h->GetCol(),h->GetRow(),h->GetQ());
  //}

  // sel hits
  for (auto h:event->GetTracks()[trkID]->GetHits()){
    if(!h->GetQ()) continue;
    event3D->Fill(h->GetTime(),h->GetRow(),h->GetCol(),h->GetQ());
    MMsel->Fill(h->GetCol(),h->GetRow(),h->GetQ());
    //MM->Fill(h->GetCol(),h->GetRow(),h->GetQ());
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
        MM->Fill(x, y, max);
    }
  }

  TCanvas *canv = new TCanvas("canv", "canv", 0., 0., 1400., 600.);
  canv->Divide(3,1);
  canv->cd(1);
  MM->Draw("COLZ");
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
  canv->WaitPrimitive();
  delete htemp;
  delete canv;

  delete MM;
  delete MMsel;
  delete event3D;
}

//******************************************************************************
std::vector<THit*> AnalysisBase::GetRobustPadsInColumn(std::vector<THit*> col) {
//******************************************************************************
  std::vector<THit*> result;
  // sort in charge decreasing order
  sort(col.begin(), col.end(), [](THit* hit1, THit* hit2){return hit1->GetQ() > hit2->GetQ();});

  for (uint i = 0; i < col.size(); ++i) {
    auto pad    = col[i];
    auto q      = pad->GetQ();
    if (!q)
      continue;

    // not more then 3 pads
    // if (i > 1)
    //   continue;

    // // WF with negative dt
    // if (pad->GetTime() - col[0]->GetTime() < -1)
    //   continue;

    // // avoid "suspisious" WF with small time difference in the 3rd pad
    // if (i > 1 && pad->GetTime() - col[0]->GetTime() < 5)
    //   continue;

    result.push_back(pad);

    // auto it_y   = pad->GetRow(_invert);
    // auto center_pad_y = geom::GetYpos(it_y, _invert);
  }

  return result;
}

//******************************************************************************
std::vector<TCluster*> AnalysisBase::GetRobustCols(std::vector<TCluster*> tr) {
//******************************************************************************
  std::vector<TCluster*> result;
  // sort clusters in increasing order
  sort(tr.begin(), tr.end(), [](TCluster* cl1,
                                TCluster* cl){
                                  return  cl1->GetCharge() < cl->GetCharge();});

  // trancation cut
  auto frac = 1.00;
  Int_t i_max = round(frac * tr.size());
  for (auto i = 0; i < i_max; ++i) {
    result.push_back(tr[i]);
  }

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
       [&](TCluster* cl1, TCluster* cl2){
          return cl1->GetX() < cl2->GetX();
        });
  return result;
}

//******************************************************************************
std::vector<TCluster*> AnalysisBase::ClusterTrack(const TTrack* tr,
                                                  int (Clustering::*f)(int, int),
                                                  Clustering& cl) {
//******************************************************************************
  std::vector<TCluster*> cluster_v;
  for (auto col:tr->GetCols(_invert)) {
    // skip first and last column
    auto it = col[0]->GetCol(_invert);
    if (it == 0 || it == geom::GetMaxColumn(_invert)-1)
      continue;
    for (auto pad:col) {
      auto col_id = pad->GetCol(_invert);
      auto row_id = pad->GetRow(_invert);

      auto cons = (cl.*f)(row_id, col_id);
      // skip first and last row
      if (row_id == 0 || row_id == geom::GetMaxRow(_invert)-1)
        continue;

      // search if the diagonal is already considered
      std::vector<TCluster*>::iterator it;
      for (it = cluster_v.begin(); it < cluster_v.end(); ++it) {
        if (!((*(*it))[0])) {
          continue;
        }

        auto cluster_col = (*(*it))[0]->GetCol(_invert);
        auto cluster_row = (*(*it))[0]->GetRow(_invert);
        if ((cl.*f)(cluster_row, cluster_col) == cons) {
          (*it)->AddHit(pad);
          (*it)->AddCharge(pad->GetQ());
          /** update X position */
          auto x_pad = geom::GetXposPad(pad, _invert, cl.angle);
          auto mult  = (*it)->GetSize();
          auto x_new = ((*it)->GetX() * (mult - 1) + x_pad) / mult;
          (*it)->SetX(x_new);
          /** */

          break;
        }
      } // loop over track clusters
      // add new cluster
      if (it == cluster_v.end()) {
        TCluster* first_cluster = new TCluster(pad);
        first_cluster->SetX(geom::GetXposPad(pad, _invert, cl.angle));
        first_cluster->SetCharge(pad->GetQ());
        cluster_v.push_back(first_cluster);
      }
    } // over pads
  } // over cols diagonalise track

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
  char *homePath(getenv("SOFTDIR"));

  if (getenv("SOFTDIR") == NULL) {
    std::cerr << "SOFTDIR varaible is not specified!" << std::endl;
    std::cerr << "Consider sourcing setup.sh" << std::endl;
    return false;
  }
  if (_param_file_name == ""){
    _param_file_name = std::string(homePath) + "/params/default.ini";
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
      auto delimiterPos = line.find("=");
      auto name = line.substr(0, delimiterPos);
      auto value = line.substr(delimiterPos + 1);
      // std::cout << name << " " << value << '\n';
      if (name == "cluster") {
        if (value == "column") {
          _clustering = CL_col;
        } else if (value == "diag") {
          _clustering = CL_diag;
          std::cout << "Diagonal cluster is used" << std::endl;
        } else if (value == "2by1") {
          _clustering = CL_2by1;
        } else if (value == "3by1") {
          _clustering = CL_3by1;
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
        } else {
          std::cout << "PRF is fit with 4th degree polinom" << std::endl;
        }
      } else if (name == "track_shape") {
        if (value == "parabola") {
          _do_para_fit = true;
          std::cout << "Parabola track fit is used" << std::endl;
        } else if (value == "linear") {
          _do_linear_fit = true;
          std::cout << "Linear track fit is used" << std::endl;
        } else {
          std::cout << "Arc track fit is used" << std::endl;
        }
      } else if (name == "max_mult") {
        _max_mult = TString(value).Atoi();
      } else if (name == "cut_gap") {
        if (value == "0") {
          _cut_gap = false;
        }
      } else if (name == "cluster_min") {
        _min_clusters = TString(value).Atoi();
      } else if (name == "max_phi") {
        _max_phi = TString(value).Atof();
      } else if (name == "max_theta") {
        _max_theta = TString(value).Atof();
//switch to WF storage
      } else if (name == "to_store_wf") {
        if (value == "0") {
          _to_store_wf = false;
        }
      }
    }
  } else {
    return false;
  }
  std::cout << "*****************************************" << std::endl;
  return true;
}

//******************************************************************************
void AnalysisBase::help(const std::string& name) {
//******************************************************************************
  std::cout << name << " usage\n" << std::endl;
  std::cout << "   -i <input_file>      : input file name with a path" << std::endl;
  std::cout << "   -o <output_path>     : output files path" << std::endl;
  std::cout << std::endl;
  std::cout << "   --start     <i>      :start from event i" << std::endl;
  std::cout << "   --end       <i>      :end with event i" << std::endl;
  std::cout << "   -t <interation>      : iteration number" << std::endl;
  std::cout << std::endl;
  std::cout << "   --param, p  <file>   : parameter file to use" << std::endl;
  std::cout << "   --prev      <file>   : file from previous iteration" << std::endl;
  std::cout << "   -b                   : run in batch mode" << std::endl;
  std::cout << "   -v <verbose_level>   : verbosity level" << std::endl;
  std::cout << "   -d                   : test mode. run over first 30 events" << std::endl;
  std::cout << "   -h                   : print ROOT help" << std::endl;
  std::cout << "   -m                   : print " << name << " help" << std::endl;
  exit(1);
}

//******************************************************************************
void AnalysisBase::process_mem_usage(double& vm_usage, double& resident_set) {
//******************************************************************************
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

//******************************************************************************
void AnalysisBase::CL_progress_dump(int eventID, int N_events) {
//******************************************************************************
  double real, virt;
  process_mem_usage(virt, real);
  double CPUtime  = _sw_event->CpuTime();
  double REALtime = _sw_event->RealTime();
  int m = 0;
  int s = 0;
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
  // std::cout << "\t Memory  " <<  real << "\t" << virt;
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
