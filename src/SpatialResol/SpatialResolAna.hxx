#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
//#include "CrossingReconstruction.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"
#include "TrackFitter.hxx"

const int Nclusters = 70;

/// Spatial resolution analysis
class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna();
  ~SpatialResolAna() override = default;

  /// Initialise histoes, input files, selections
  bool Initialize(int argc, char** argv) override;

  /// Read CLI
  bool ReadCLI(int argc, char **argv) override;

  /// Process the selection output called Event
  bool ProcessEvent(const TEvent* event) override;

  /// Draw the histograms of interest
  TCanvas* DrawSelectionCan(const TRawEvent* event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput() override;

  /// Profile PRF with peak and RMS
  bool ProfilePRF(const TH2F* PRF_h, TGraphErrors* gr);
  /// Profile PRF along X axis
  bool ProfilePRF_X(const TH2F* PRF_h, TGraphErrors* gr, TH1F* errors);

  /// Initialise PRF with expected params
  TF1* InitializePRF(TString name, bool shift=false);

  /// Get mean and FWHM for the histo
  Double_t GetFWHM(const TH1F* h, Double_t& mean);
  Double_t GetFWHM(const TH1F* h);

  /// Draw the histograms of interest
  bool Draw();

  /// verbosity levels
  enum verbosity_SR {
    v_analysis_steps = v_base_last + 1,
    v_fit_details,
    v_residuals,
    v_prf
  };

 protected:
  /// Previous iteration output to extract PRF
  TFile*  _prev_iter_file;

  /// Name of the file from previous iteration
  TString _prev_iter_name;
  /// iteration number. Starting from 0
  Int_t   _iteration;

  /// Whether to apply correction of spatial resolution (take geometrical mean)
  bool _correction;

  /// output tree
  TTree*  _tree;
  /// Oputput tree vars
  // Event vars
  /// event number
  Int_t   _ev;
  /// angle in MM plane
  Float_t _angle_xy;
  /// angle w.r.t. MM
  Float_t _angle_yz;

  /// Number of robust clusters in event
  Int_t _rob_clusters;

  /// track fit quality Chi2/NDF
  Float_t _quality;
  /// momentum
  Float_t _mom;
  /// angle
  Float_t _sin_alpha;
  /// offset
  Float_t _offset;

  /// Cluster vars
  /// Position of the cluster
  Float_t _clust_pos[Nclusters];
  /// X position of the cluster
  Float_t _x[Nclusters];
  /// Y Position of the cluster
  // Float_t _cluster_av[Nclusters];
  /// X position of the avareged cluster
  Float_t _x_av[Nclusters];
  /// Position of the track
  Float_t _track_pos[Nclusters];
  /// Residuals (X_track-X_cluster)
  Float_t _residual[Nclusters];
  /// Residuals (X_track-X_cluster) w/o the given cluster in the fit
  Float_t _residual_corr[Nclusters];
  /// charge in the cluster
  Int_t   _charge[Nclusters];
  /// multiplicity of the cluster
  Int_t   _multiplicity[Nclusters];

  /// Track analytical fit function for plotting
  TF1* _track_fit_func;

  /// Pad vars
  /// X_track - X_pad --> X axis of the PRF
  Float_t _dx[Nclusters][10];
  /// Fraction of charge Q_pad / Q_cluster --> Y axis of PRF
  Float_t _qfrac[Nclusters][10];
  /// time of the pad
  Int_t   _time[Nclusters][10];

  /** variables for the dEdx analysis */
  /// dE/dx
  Float_t _dEdx;
  TH1F* _hdEdx;
  Int_t _pad_charge[Nclusters][10];
  Int_t _pad_time[Nclusters][10];

  /// pad X position == column
  Int_t _pad_x[Nclusters][10];
  /// pad Y position === row
  Int_t _pad_y[Nclusters][10];

  /// Waveform width
  Int_t _wf_width[Nclusters][10];
  /// Waveform width at half maximum
  Int_t _wf_fwhm[Nclusters][10];


  //Int_t _pad_wf_t[Nclusters][10][520];
  Int_t _pad_wf_q[Nclusters][10][520];

  /** Histograms **/
  /** Pad response function block **/
  /// PRF function from the previous step. Used for Chi2 fit
  TF1*  _prf_function;
  TF1**  _prf_function_arr;
  /// PRF histoes
  TH2F* _prf_histo;
  // PRF profiling graphs
  TGraphErrors* _prf_graph;

  // TF1*    _prf_function_2pad;
  // TF1*    _prf_function_3pad;
  // TF1*    _prf_function_4pad;

  TH2F* _prf_histo_2pad;
  TH2F* _prf_histo_3pad;
  TH2F* _prf_histo_4pad;

  TGraphErrors* _prf_graph_2pad;
  TGraphErrors* _prf_graph_3pad;
  TGraphErrors* _prf_graph_4pad;

  /// Pad response function in time
  TH2F* _prf_time;

  /// analytical PRF time function
  TF1* _prf_time_func;

  /// uncertainties of the time profile
  TGraphErrors* _prf_time_error;

  TH1F* _prf_time_e;

  /// WARNING TEMP
  Float_t _fit_up[Nclusters];
  Float_t _fit_bt[Nclusters];

  /// Fitter class for the track and cluster fitting
  TrackFitCern* _fitter;

  /** Switchers **/
  /// Whether to use full track fitting
  bool _do_full_track_fit;

  /// Whether fit all the pads separatly
  bool _do_separate_pad_fit;

  /// Whether to fit residuals with Gaussian
  bool _gaussian_residuals;

  /// Whether to assign uncertainty to charge
  bool _charge_uncertainty;

  TH1F*   _mom_reco;
  TH1F*   _pos_reco;
  TH1F*   _ang_reco;
  TH1F*   _qulity_ratio;
  TH1F*   _chi2_ratio;

  /// Chi2 function of the track fit
  TH1F* _chi2_track;

  /// How many columns are used for fit
  TH1F* _cols_used;

  /// Residuals X_track - X_fit histoes
  TH1F* _resol_total;

  TH1F* _resol_col_hist[Nclusters];
  TH1F* _resol_col_hist_except[Nclusters];

  TH1F* _resol_col_hist_2pad[Nclusters];
  TH1F* _resol_col_hist_2pad_except[Nclusters];

  TH1F* _resol_col_hist_3pad[Nclusters];
  TH1F* _resol_col_hist_3pad_except[Nclusters];

  TGraphErrors* _residual_mean;
  TGraphErrors* _residual_sigma;

  TGraphErrors* _residual_sigma_unbiased;
  TGraphErrors* _residual_sigma_biased;

  TH2F* _prf_histo_col[Nclusters];

  /// separate pad fit study
  TH1F* _fit_quality_plots[3][Nclusters];
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
  TH1F* _resol_col_x_scan[Nclusters][x_scan_bin];
  TH1F* _mult_x_scan[Nclusters][x_scan_bin];
//  TH1F* _x_pads = new TH1F("padX", "", 4, -0.03, 0.01);

  TH1F* _resol_col_x_scan_lim_mult[Nclusters][x_scan_bin];

  TH2F* _prf_histo_xscan[4];
  TGraphErrors* _prf_graph_xscan[4];

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

  /// space limit for the PRF usage
  /// pads that are far away are supposed to be unreliable
  const float fit_bound_left  = -0.025;
  const float fit_bound_right =  0.025;
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
