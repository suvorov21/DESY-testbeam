#include <iostream>  // stream
#include <unistd.h>  // getopt on Mac
#include <fstream>   // read file lists
#include <algorithm>

#include "TROOT.h"

#include "AnalysisBase.hxx"

AnalysisBase::AnalysisBase(int argc, char** argv) {
  // default values
  _verbose    = 1;
  _batch      = false;
  _test_mode  = false;

  _file_in_name           = "";
  _file_out_name          = "";
  _event_list_file_name   = "";

  _file_in    = NULL;
  _file_out   = NULL;

  _chain      = NULL;
  _reconstruction  = NULL;

  // read CLI
  for (;;) {
    int c = getopt(argc, argv, "i:o:bv:dmt:");
    if (c < 0) break;
    switch (c) {
      case 'i' : _file_in_name     = optarg;       break;
      case 'o' : _file_out_name    = optarg;       break;
      case 'b' : _batch            = true;         break;
      case 'v' : _verbose          = atoi(optarg); break;
      case 'd' : _test_mode        = true;         break;
      case 'm' : help(argv[0]);                    break;
      case 't' : _iteration        = atoi(optarg); break;
      case '?' : help(argv[0]);
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
  std::cout << "Initializing analysis base...............";
  // read and chain input files
  _chain = new TChain("tree");

  if (_file_in_name.Contains(".root")) {
    _chain->AddFile(_file_in_name);
  } else {
    std::ifstream fList(_file_in_name.Data());
    if (!fList.good()) {
      std::cerr << "Can not read input " << _file_in_name << std::endl;
      exit(1);
    }
    while (fList.good()) {
      std::string filename;
      getline(fList, filename);
      if (fList.eof()) break;
      _chain->AddFile(filename.c_str());
    }
  }

  _chain->SetBranchAddress("PadAmpl", _padAmpl);

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

  if (_event_list_file_name == "") {
    Int_t N_events = _chain->GetEntries();

    for (auto i = 0; i < N_events; ++i)
      _EventList.push_back(i);
  } else {
    // FIXIT read the event list file
  }

  // Open the output file
  _file_out = new TFile(_file_out_name.Data(), "NEW");
  if (!_file_out->IsOpen()) {
    std::cerr << "ERROR. AnalysisBase::Initialize()" << std::endl;
    std::cerr << "File already exists or directory is not writable" << std::endl;
    std::cerr << "To prevent overwriting of the previous result the program will exit" << std::endl;
    exit(1);
  }

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
    N_events = std::min(static_cast<Int_t>(EventList.size()), 1000);

  if (_verbose == 1) {
    std::cout << "Processing" << std::endl;
    std::cout << "[                               ]   Nevents = " << N_events << "\r[";
  }

  for (auto eventID = 0; eventID < N_events; ++eventID) {
    if (_verbose > 1)
      std::cout << "Event " << eventID << std::endl;

    if (_verbose == 1 && (eventID%(N_events/30)) == 0)
      std::cout << "." << std::flush;

    _chain->GetEntry(EventList[eventID]);

    Event event;
    event.trackNum = 0;

    if (!_reconstruction->SelectEvent(_padAmpl, event))
      continue;

    ProcessEvent(event);
  }

  if (_verbose == 1)
    std::cout << "]" << std::endl;

  return true;
}

bool AnalysisBase::ProcessEvent(const Event event) {
  (void)event;
  std::cerr << "EROOR. AnalysisBase::ProcessEvent(). Event processing should be defined in your analysis" << std::endl;
  exit(1);
  return true;
}

bool AnalysisBase::WriteOutput() {
  // WARNING add error

  if (!_file_out->IsOpen())
    return false;

  std::cout << "Writing standard output..................";

  _file_out->cd();

  auto size = static_cast<int>(_output_vector.size());
  for (auto i = 0; i < size; ++i)
    _output_vector[i]->Write();

  _file_out->Close();

  std::cout << "done     " << "Write  " << size << " objects" << std::endl;

  return true;
}

void AnalysisBase::help(const std::string name) {
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

