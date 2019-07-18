#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TApplication.h"

#include "ReconstructionBase.hxx"
#include "SetT2KStyle.hxx"

/// Main analysis template
class AnalysisBase {
 public:
  AnalysisBase(int argc, char** argv);
  virtual ~AnalysisBase() {;}

  /// Initialise histoes, input files, selections
  virtual bool Initialize();
  /// Loop over TChain entries. Can use pre-defined event list
  virtual bool Loop(std::vector<Int_t> EventList);
  /// Process the selection output called Event
  virtual bool ProcessEvent(const Event event);
  /// Write output files (histos, trees)
  virtual bool WriteOutput();

  /// Print usage
  void help(const std::string name);

  std::vector<Int_t> GetEventList() {return _EventList;}

 protected:
  TString _file_in_name;
  TString _file_out_name;

  TString _event_list_file_name;
  std::vector<Int_t> _EventList;

  TFile* _file_in;
  TFile* _file_out;

  TChain* _chain;

  /// what we read from input
  Int_t _padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];

  /// outout vector to put in the file
  std::vector<TObject*> _output_vector;

  /// Selection. You can use plenty in the analysis.
  /** At least one should be defines */
  ReconstructionBase* _selection;

  /// T2K plotting style
  TStyle* _t2kstyle;

  /// DEBUG vars
  Int_t _verbose;
  bool _batch;
  bool _test_mode;

  TApplication* _app;
};


#endif  // SRC_BASE_ANALYSISBASE_HXX_
