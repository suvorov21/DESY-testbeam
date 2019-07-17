#ifndef SRC_DEDX_DEDXANA_HXX_
#define SRC_DEDX_DEDXANA_HXX_

#include "AnalysisBase.hxx"
#include "DBSCANSelection.hxx"

/// Spatial resolution analysis
class dEdxAna: public AnalysisBase {
 public:
  dEdxAna(int argc, char** argv);
  virtual ~dEdxAna() {;}

  TH1F* _hdEdx;
  int   _selEvents;

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const Event event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
