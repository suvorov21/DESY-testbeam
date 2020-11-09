#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
//#include "CrossingReconstruction.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"
#include "TrackFitter.hxx"

/// Spatial resolution analysis
class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna(int argc, char** argv);
  virtual ~SpatialResolAna() {;}

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const TEvent* event);

  /// Draw the histograms of interest
  TCanvas* DrawSelectionCan(const TEvent* event, int trkID);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();

  bool ProfilePRF(const TH2F* _PRF_h, TGraphErrors* gr);

  TF1* InitializePRF(const TString name);

  Double_t GetFWHM(const TH1F* h, Double_t& mean);

 private:
  /// Previous iteration output to extract PRF
  TFile*  _Prev_iter_file;

  /// output tree
  TTree* _tree;
  Float_t _angle_xy;
  Float_t _angle_yz;
  Float_t _residual[geom::nPadx];
  Int_t   _charge[geom::nPadx];
  Int_t   _multiplicity[geom::nPadx];
  Float_t _dx[geom::nPadx][10];
  Float_t _qfrac[geom::nPadx][10];

  /// PRF function from the previous step. Used for Chi2 fit
  TF1*    _PRF_function;
  /// PRF histoes
  TH2F* _PRF_histo;
  // PRF profiling graphs
  TGraphErrors* _PRF_graph;

  TF1*    _PRF_function_2pad;
  TF1*    _PRF_function_3pad;
  TF1*    _PRF_function_4pad;

  TH2F* _PRF_histo_2pad;
  TH2F* _PRF_histo_3pad;
  TH2F* _PRF_histo_4pad;

  TGraphErrors* _PRF_graph_2pad;
  TGraphErrors* _PRF_graph_3pad;
  TGraphErrors* _PRF_graph_4pad;

  /// Fitter class for the track and cluster fitting
  TrackFitter* _fitter;

  /// Whether to use arc function for track fitting
  bool _do_arc_fit;
  /// Whether to use full track fitting
  bool _do_full_track_fit;

  /// Whether fit all the pads separatly
  bool _do_separate_pad_fit;

  /// Whether to apply correction of spatial resolution (take geometrical mean)
  bool _correction;

  /// Whether to fit residuals with Gaussian
  bool _gaussian_residuals;

  /// Whether to assign uncertainty to charge
  bool _charge_uncertainty;

  /// Whether to use Gaussian lorentzian PRf fit over polynomial
  bool _gaus_lorentz_PRF;

  /// iteration number. Starting from 0
  Int_t   _iteration;

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
  TH1F* _resol_total;

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

  TH2F* _PRF_histo_col[geom::nPadx];

  /// separate pad fit study
  TH1F* _Fit_quality_plots[3][geom::nPadx];
  TAxis* _prf_scale_axis;

  /// errors vs the PRF value
  static const int prf_error_bins = 10;
  Double_t prf_error_bins_arr[prf_error_bins] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01};
  TAxis* _prf_error_axis = new TAxis(prf_error_bins-1, prf_error_bins_arr);
  TH1F* _uncertainty_prf_bins[prf_error_bins-1];
  TGraphErrors* _uncertainty_vs_prf_gr;
  TGraphErrors* _uncertainty_vs_prf_gr_prev;
  TH1F*         _uncertainty_vs_prf_histo;


  /// x scan data
  static const int x_scan_bin = 50;
  const float x_scan_min = -0.035;
  const float x_scan_max = 0.015;
  TAxis* _x_scan_axis;
  TH1F* _resol_col_x_scan[geom::nPadx][x_scan_bin];
  TH1F* _mult_x_scan[geom::nPadx][x_scan_bin];
  TH1F* _x_pads = new TH1F("padX", "", 4, -0.03, 0.01);

  TH1F* _resol_col_x_scan_lim_mult[geom::nPadx][x_scan_bin];

  TH2F* _PRF_histo_xscan[4];
  TGraphErrors* _PRF_graph_xscan[4];

  /// Average uncertainty from the previous iteration
  Float_t _uncertainty;

  /// vector of events IDs that passed the Reco and selection
  std::vector<Int_t> _passed_events;

  // [units are meters]
  const float prf_min     = -0.027;
  const float prf_max     = 0.027;
  const int   prf_bin     = 180;

  const float resol_min   = -0.004;
  const float resol_max   = 0.004;
  const int   resol_bin   = 200.;

  const float fit_bound_left  = -0.025;
  const float fit_bound_right =  0.025;
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
