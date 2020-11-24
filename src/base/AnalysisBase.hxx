#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

/** @cond */
#include <numeric>

#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <TNtuple.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
#include "TStopwatch.h"
#include "TGraphErrors.h"
/** @endcond */

#include "ReconstructionBase.hxx"
#include "SetT2KStyle.hxx"

class TCluster;

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

  virtual void CL_progress_dump(int eventID, int Nevents);

  std::vector<THit*> GetRobustPadsInColumn(std::vector<THit*> col);
  std::vector<TCluster*> GetRobustCols(std::vector<TCluster*> tr);

  /// Print usage
  void help(const std::string& name);

  void process_mem_usage(double& vm_usage, double& resident_set);

  void SetEventList(const std::vector<Int_t>& var) {_EventList.clear(); _EventList = var;}
  std::vector<Int_t> GetEventList() const {return _EventList;}

  /// Clusterization
  std::vector<TCluster*> DiagonolizeTrack(const TTrack* tr);

  std::vector<TCluster*> ColonizeTrack(const TTrack* tr);

  AnalysisBase(const AnalysisBase& ana){(void)ana;
    std::cerr << "Copy constructor is depricated" << std::endl; exit(1);}
  bool operator==(const AnalysisBase* ana){(void)ana;
    std::cerr << "Comparison is depricated" << std::endl; exit(1);}

 protected:

  TString _file_in_name;
  TString _file_out_name;

  TString _event_list_file_name;
  std::vector<Int_t> _EventList;
  bool    _store_event;
  int     _start_ID;
  int     _end_ID;
  int     _selected;

  TEvent* _event;
  bool    _store_event_tree;
  TFile*  _event_file;
  TTree*  _event_tree;
  bool    _work_with_event_file;

  TFile* _file_in;
  TFile* _file_out;

  TChain* _chain;

  /// what we read from input
  bool _saclay_cosmics;
  Int_t _padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];
  Int_t _padAmpl_saclay[geom::nPadx][geom::nPady][geom::Nsamples_saclay];

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

  /// Whether to invert track analysis logic
  /// E.g. analyse cosmic tracks
  bool _invert;

  TApplication* _app;
  TStopwatch* _sw_event;

  TStopwatch* _sw_partial[5];

  /// Use CERN data
  bool _useCern;
  TTree* _tgeom;

  std::vector<short>          *_listOfChannels;
  std::vector<std::vector<short> > *_listOfSamples;

  std::vector<int> *_iPad;
  std::vector<int> *_jPad;
};


#endif  // SRC_BASE_ANALYSISBASE_HXX_
