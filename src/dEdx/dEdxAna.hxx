#ifndef SRC_DEDX_DEDXANA_HXX_
#define SRC_DEDX_DEDXANA_HXX_

#include "AnalysisBase.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"

const double alpha = 0.625;


/// Spatial resolution analysis
class dEdxAna: public AnalysisBase {
 public:
  dEdxAna(int argc, char** argv);
  virtual ~dEdxAna() {;}

  /**
   * Trees declaration
   */
  /// Tree with entry per event
  TTree *outtree;
  Int_t _ev;
  Float_t _dEdx;
  Float_t _angle_xy;
  Float_t _angle_yz;
  Int_t _npoints;

  Int_t _multiplicity[geom::nPadx];
  Int_t _multiplicity_robust[geom::nPadx];
  Int_t _charge[geom::nPadx];
  Int_t _maxcharge_time[geom::nPadx];
  Float_t _maxcharge_frac[geom::nPadx];

  Int_t _pad_charge[10][geom::nPadx];
  Int_t _pad_time[10][geom::nPadx];

  Int_t _pad_x[10][geom::nPadx];
  Int_t _pad_y[10][geom::nPadx];

  Int_t _wf_width[10][geom::nPadx];
  Int_t _wf_fwhm[10][geom::nPadx];

  /**
   * Histogramms and graphs declaration
   */

  TH1F* _hdEdx;
  TH1F* _hTime;
  TH1F* _mult;
  Int_t   _selEvents;

  TGraphErrors* _mult_graph;
  TGraphErrors* _mult_graph_err;

  /// maximum charge in pad in event
  TH1F* _max_charge_pad;

  /// cluster charge before truncation
  TH1F* _un_trunk_cluster;

  /// Max charge pos and time
  TH1F* _max_charge_pos;
  TH1F* _max_charge_time;

  /// Angular histo
  TH2F* _angle;

  /// Time difference between leading and subleading pads
  TH1F* _delta_t_fst;
  TH1F* _delta_t_scd;

  /// Time difference wrt angle in YZ
  TH2F* _delta_t_angle;

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
