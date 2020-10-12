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
  TGraphErrors* _mult_graph_err;

  /// maximum charge in pad in event
  TH1F* _max_charge_pad;

  /// cluster charge before truncation
  TH1F* _un_trunk_cluster;
  /// leading pad charge
  TH1F* _fst_pad_charge;
  TH1F* _scd_pad_charge;
  TH1F* _trd_pad_charge;
  TH1F* _fth_pad_charge;

  /// Max charge pos and time
  TH1F* _max_charge_pos;
  TH1F* _max_charge_time;

  /// Angular histo
  TH2F* _angle;

  ///
  std::vector<TH1F*> _charge_per_mult;

  /// HIT MAPS
  TH2F* _XZ_leading;
  TH1F* _XZ_bias;

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const TEvent *event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();

  ///
  bool DrawCharge();
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
