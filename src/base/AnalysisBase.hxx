#ifndef ANALYSISBASE_H
#define ANALYSISBASE_H

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TApplication.h"

#include "../../utils/Geom.hxx"
#include "SelectionBase.hxx"

class AnalysisBase {
  AnalysisBase(int argc, char** argv);
  virtual ~AnalysisBase() {}

public:
  // Initialise histoes, input files, selections
  bool Initialise();
  // loop over TChain
  bool Loop();
  // Process the selection output called Event
  bool ProcessEvent(const Event event);
  // write output files (histos, trees)
  bool WriteOutput();

  // print usage
  void help(std::string name);

protected:
  TString _file_in_name;
  TString _file_out_name;

  TFile* _file_in;
  TFile* _file_out;

  TChain* _chain;

  // what we read from input
  Int_t _padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];

  // Selection. You can use plenty in the analysis.
  // At least one should be defines
  SelectionBase* _selection;

  // DEBUG vars
  Int_t _verbose;
  bool _batch;
  bool _test_mode;

  TApplication* _app;
};


#endif