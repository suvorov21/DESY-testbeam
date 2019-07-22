#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

/** @cond */
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <TNtuple.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
/** @endcond */

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
  virtual bool ProcessEvent(const TEvent* event);
  /// Write output files (histos, trees)
  virtual bool WriteOutput();
  virtual void DrawSelection(const TEvent *event, int trackID);

  /// Print usage
  void help(const std::string name);

  void process_mem_usage(double& vm_usage, double& resident_set);

  void SetEventList(std::vector<Int_t> var) {_EventList = var;}
  std::vector<Int_t> GetEventList() const {return _EventList;}

 protected:
  /// iteration number. Starting from 0
  // TODO remove it out to particular analysis
  // Don't know how to parse CLI in different classes
  Int_t   _iteration;

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
  ReconstructionBase* _reconstruction;

  /// T2K plotting style
  TStyle* _t2kstyle;

  /// DEBUG vars
  Int_t _verbose;
  bool _batch;
  bool _test_mode;
  bool _overwrite;

  TApplication* _app;
};


#endif  // SRC_BASE_ANALYSISBASE_HXX_
