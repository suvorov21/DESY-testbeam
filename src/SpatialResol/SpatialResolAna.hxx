#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
//#include "CrossingReconstruction.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"

/// Spatial resolution analysis
class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna(int argc, char** argv);
  virtual ~SpatialResolAna() {;}

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const TEvent* event);

  /// Fit the whole track with CERN method
  TF1* GetTrackFitCERN(const double* track_pos, const int* mult, const int miss_id = -1);
  /// Fit the whole track with ILC method
  TF1* GetTrackFitILC(const TTrack* track, const double pos, const int miss_id = -1);

  /// Extract cluster position with CERN method
  double GetClusterPosCERN(const std::vector<THit*>& col, const int cluster, const double pos);
  /// Extract cluster position with ILC method
  double GetClusterPosILC(const std::vector<THit*>& col, const double pos);

  /// Whether to miss the column in the fitter
  bool MissColumn(int it_x);
  /// Whether the cluster is good for fitting
  bool UseCluster(const std::vector<THit*>& col);

  /// Get the number of used columns
  Int_t GetMaxColumn();

  /// Draw the histograms of interest
  TCanvas* DrawSelectionCan(const TEvent* event, int trkID);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();
 private:
  /// Previous iteration output to extract PRF
  TFile*  _Prev_iter_file;
  /// PRF function from the previous step. Used for Chi2 fit
  TF1*    _PRF_function;

  /// Whether to use arc function for track fitting
  bool    _do_arc_fit;
  /// Whether to use full track fitting
  bool    _do_full_track_fit;

  /// Whether to apply correction of spatial resolution (take geometrical mean)
  bool    _correction;

  /// Whether to fit residuals with Gaussian
  bool _gaussian_residuals;

  /// Whether to assign uncertainty to charge
  bool _charge_uncertainty;

  /// iteration number. Starting from 0
  Int_t   _iteration;

  /// Whether to invert track analysis logic
  /// E.g. analyse cosmic tracks
  bool _invert;

  /// Fitting function for track going up
  TF1*    _circle_function_up;
  /// Fitting function for track going down
  TF1*    _circle_function_dn;

  TH1F*   _mom_reco;
  TH1F*   _pos_reco;
  TH1F*   _ang_reco;
  TH1F*   _qulity_ratio;
  TH1F*   _chi2_ratio;

  /// Chi2 function of the track fit
  TH1F* _Chi2_track;

  /// How many columns are used for fit
  TH1F* _Cols_used;

  /// Residuals X_track - X_fit histoes
  TH1F* _resol_col_hist[geom::nPadx];
  TH1F* _resol_col_hist_except[geom::nPadx];

  TH1F* _resol_col_hist_2pad[geom::nPadx];
  TH1F* _resol_col_hist_2pad_except[geom::nPadx];

  TH1F* _resol_col_hist_3pad[geom::nPadx];
  TH1F* _resol_col_hist_3pad_except[geom::nPadx];

  TGraphErrors* _residual_mean;
  TGraphErrors* _residual_sigma;
  // TH1F* _residual_sigma_2pad;
  // TH1F* _residual_sigma_3pad;

  TGraphErrors* _residual_sigma_unbiased;
  TGraphErrors* _residual_sigma_biased;

  /// PRF histoes
  TH2F* _PRF_histo;
  TH2F* _PRF_histo_2pad;
  TH2F* _PRF_histo_3pad;
  TH2F* _PRF_histo_4pad;

  TH2F* _PRF_histo_col[geom::nPadx];

  static const int x_scan_bin = 50;
  const float x_scan_min = -0.035;
  const float x_scan_max = 0.015;
  TAxis* _x_scan_axis;
  TH1F* _resol_col_x_scan[geom::nPadx][x_scan_bin];
  TH1F* _mult_x_scan[geom::nPadx][x_scan_bin];

  TH1F* _resol_col_x_scan_lim_mult[geom::nPadx][x_scan_bin];

  /// Average uncertainty from the previous iteration
  Float_t _uncertainty;

  // PRF profiling graphs
  TGraphErrors* _PRF_graph;

  /// vector of events IDs that passed the Reco and selection
  std::vector<Int_t> _passed_events;

  // [units are meters]
  const float prf_min     = -0.027;
  const float prf_max     = 0.027;
  const int   prf_bin     = 180;

  const float resol_min   = -0.004;
  const float resol_max   = 0.004;
  const int   resol_bin   = 200.;

  const float fit_bound_left  = -0.015;
  const float fit_bound_right =  0.015;

  const float default_error   = 0.001;
  const float one_pad_error   = 0.002;
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
