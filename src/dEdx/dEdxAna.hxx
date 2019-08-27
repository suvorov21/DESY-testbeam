#ifndef SRC_DEDX_DEDXANA_HXX_
#define SRC_DEDX_DEDXANA_HXX_

#include "AnalysisBase.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"


/// Spatial resolution analysis
class dEdxAna: public AnalysisBase {
 public:
  dEdxAna(int argc, char** argv);
  virtual ~dEdxAna() {;}

  TH1F* _hdEdx;
  TH1F* _hTime;
  TH1F* _mult;
  int   _selEvents;

  TH1F* _mult_col[geom::nPadx];
  TGraphErrors* _mult_graph;

  TH1F* _max_charge_pad;

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const TEvent *event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
