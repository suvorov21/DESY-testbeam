#include "AnalysisBase.hxx"

#include <iostream> // stream
#include <unistd.h> // getopt on Mac
#include <fstream>  // read file lists

#include "TROOT.h"

#include "../../utils/SetT2KStyle.hxx"

AnalysisBase::AnalysisBase(int argc, char** argv) {
  // default values
  _verbose    = 1;
  _batch      = false;
  // read CLI
  for (;;) {
    int c = getopt(argc, argv, "i:o:bvdm");
    if (c < 0) break;
    switch (c) {
      case 'i' : _file_in_name     = optarg;       break;
      case 'o' : _file_out_name    = optarg;       break;
      case 'b' : _batch            = true;         break;
      case 'v' : _verbose          = atoi(optarg); break;
      case 'd' : _test_mode        = true;         break;
      case 'm' : help(argv[0]);                    break;
      case '?' : help(argv[0]);
    }
  }

  if (!_batch)
    _app = new TApplication("app", &argc, argv);
}

bool AnalysisBase::Initialise() {
  // read and chain input files
  auto chain = new TChain("tree");

  if (_file_in_name.Contains(".root")) {
    std::cout << "adding filename" <<" " << _file_in_name << std::endl;
    chain->AddFile(_file_in_name);
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
      chain->AddFile(filename.c_str());
    }
  }

  chain->SetBranchAddress("PadAmpl", _padAmpl);

  // setup the T2K style
  auto T2KstyleIndex = 2;
  // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
  auto localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper

  auto t2kstyle = SetT2KStyle(T2KstyleIndex, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());
  gROOT->ForceStyle();

  // Initilise selection
  _selection = new SelectionBase();

  // Initialise histoes
  // * do it in your analysis *

  return true;
}

bool AnalysisBase::Loop() {
  Int_t N_events = _chain->GetEntries();
  if (_test_mode)
    N_events = std::min(N_events, 30);

  if (_verbose == 1) {
    std::cout << "Processing" << std::endl;
    std::cout << "[                    ]   Nevents = " << N_events << "\r[";
  }


  for (auto eventID = 0; eventID < N_events; ++eventID) {
    if (_verbose > 1)
      std::cout << "Event " << eventID << std::endl;

    if (_verbose == 1 && (eventID%(N_events/20)) == 0)
      std::cout << "." << std::flush;

    _chain->GetEntry(eventID);

    Event event;

    if (!_selection->SelectEvent(_padAmpl, event))
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

void AnalysisBase::help(std::string name) {
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

