#include <algorithm>

#include "TROOT.h"

#include "AnalysisBase.hxx"

AnalysisBase::AnalysisBase(int argc, char** argv) :
  _file_in_name(""),
  _file_out_name(""),
  _event_list_file_name(""),
  _event(NULL),
  _start_ID(-1),
  _end_ID(-1),
  _selected(0),
  _store_event_tree(false),
  _work_with_event_file(false),
  _file_in(NULL),
  _file_out(NULL),
  _chain(NULL),
  _reconstruction(NULL),
  _verbose(1),
  _batch(false),
  _test_mode(false),
  _overwrite(false),
  _invert(false),
  _app(NULL),
  _useCern(false)
{

  // TODO redefine CLI.
  // now written in an ugly way
  // prevent copy-past between the daughter-parent classes
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

  // read CLI
  for (;;) {
    int c = getopt_long(argc, argv, "i:o:bv:drhst:ca", longopts, &index);
    if (c < 0) break;
    switch (c) {
      case 0  :
        if (index == 11) _start_ID         =  atoi(optarg);
        if (index == 12) _end_ID           =  atoi(optarg);
        break;
      case 'i' : _file_in_name     = optarg;       break;
      case 'o' : _file_out_name    = optarg;       break;
      case 'b' : _batch            = true;         break;
      case 'v' : _verbose          = atoi(optarg); break;
      case 'd' : _test_mode        = true;         break;
      case 'r' :
        _overwrite        = true;
        std::cout << "Output will be overwritten" << std::endl;
        break;
      case 'h' : help(argv[0]);                    break;
      case 'a' :
        _invert           = true;
        std::cout << "Inverted geometry is used" << std::endl;
        break;
      case 's' :
        _store_event_tree = true;
        std::cout << "Tree with TEvents will be written" << std::endl;
        break;
      //case '?' : help(argv[0]);
    }
  }

  if (_file_in_name == "") {
    std::cerr << "ERROR. AnalysisBase::AnalysisBase. No input file specified" << std::endl;
    exit(1);
  }

  if (!_batch)
    _app = new TApplication("app", &argc, argv);
}

bool AnalysisBase::Initialize() {

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
    tree_name = "padData";
    std::cout << "Running over CERN data" << std::endl;
    _tgeom= (TTree*)file->Get("femGeomTree");
    //std::vector<int> *iPad(0);
    //std::vector<int> *jPad(0);
    _tgeom->SetBranchAddress("jPad", &_jPad );
    _tgeom->SetBranchAddress("iPad", &_iPad );
    _tgeom->GetEntry(0); // put into memory geometry info
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

    _chain->SetBranchAddress("PadphysChannels", &_listOfChannels );
    _chain->SetBranchAddress("PadADCvsTime"   , &_listOfSamples );

  } else if (!_work_with_event_file) {
    _chain->SetBranchAddress("PadAmpl", _padAmpl);
    // read the branch name in format "padAmple[][][]"
    TString branch_name = _chain->GetBranch("PadAmpl")->GetTitle();
    if ((branch_name.Contains("[510]") && geom::Nsamples != 510) ||
        (branch_name.Contains("[511]") && geom::Nsamples != 511)) {
      std::cerr << "ERROR. AnalysisBase::Initialize()" << std::endl;
      std::cerr << "Time binning is inconsistent." << std::endl;
      std::cerr << "File branch is " << branch_name << " while " << geom::Nsamples << " is read" << std::endl;
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

bool AnalysisBase::Loop(std::vector<Int_t> EventList) {
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

  if (_verbose == 1) {
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
    if (_verbose > 1)
      std::cout << "Event " << eventID << std::endl;

    if (_verbose == 1 && (eventID%(N_events/denimonator)) == 0)
      this->CL_progress_dump(eventID, N_events);

    _chain->GetEntry(EventList[eventID]);

    if (_useCern) {
      memset(_padAmpl, 0, geom::nPadx * geom::nPady * geom::Nsamples * (sizeof(Int_t)));
      for (uint ic=0; ic< _listOfChannels->size(); ic++){
        int chan= (*_listOfChannels)[ic];
        if ((*_jPad)[chan] >= geom::nPadx ||
              (*_iPad)[chan] >= geom::nPady)
          continue;
        for (uint it = 0; it < (*_listOfSamples)[ic].size(); it++){
          if (it >= geom::Nsamples)
            continue;
          int adc = (*_listOfSamples)[ic][it];
          _padAmpl[(*_jPad)[chan]][(*_iPad)[chan]][it] = adc;
        }
      }
    } else {
      // Subtract the pedestal
      for (auto x = 0; x < geom::nPadx; ++x) {
        for (auto y = 0; y < geom::nPady; ++y) {
          for (auto t = 0; t < geom::Nsamples; ++t) {
            if (_padAmpl[x][y][t] - 250 < 0)
              _padAmpl[x][y][t] = 0;
            else
              _padAmpl[x][y][t] = _padAmpl[x][y][t] - 250;
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
    else _event->SetID(EventList[eventID]);

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

  if (_verbose == 1)
    std::cout << std::endl;

  return true;
}

bool AnalysisBase::ProcessEvent(const TEvent* event) {
  (void)event;
  std::cerr << "EROOR. AnalysisBase::ProcessEvent(). Event processing should be defined in your analysis" << std::endl;
  exit(1);
  return true;
}

bool AnalysisBase::WriteOutput() {

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
void AnalysisBase::DrawSelection(const TEvent *event, int trkID){
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

std::vector<THit*> AnalysisBase::GetRobustPadsInColumn(std::vector<THit*> col) {
  std::vector<THit*> result;
  // sort in charge decreasing order
  sort(col.begin(), col.end(), [](THit* hit1, THit* hit2){return hit1->GetQ() > hit2->GetQ();});

  int cluster_q = 0;
  // for (auto pad:col)
  //   cluster_q += pad->GetQ();

  // if (cluster_q > 1800)
  //   return result;

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

void AnalysisBase::help(const std::string& name) {
  std::cout << name << " usage\n" << std::endl;
  std::cout << "   -i <input_file>      : input file name with a path" << std::endl;
  std::cout << "   -o <output_path>     : output files path" << std::endl;
  std::cout << "   -b                   : run in batch mode" << std::endl;
  std::cout << "   -v <verbose_lvel>    : verbosity level" << std::endl;
  std::cout << "   -d                   : test mode. run over first 30 events" << std::endl;
  std::cout << "   -h                   : print ROOT help" << std::endl;
  std::cout << "   -m                   : print " << name << " help" << std::endl;
  exit(1);
}

void AnalysisBase::process_mem_usage(double& vm_usage, double& resident_set)
{
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

void AnalysisBase::CL_progress_dump(int eventID, int N_events) {
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
