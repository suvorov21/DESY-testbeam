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

  /// Profile PRF with peak and RMS
  bool ProfilePRF(const TH2F* _PRF_h, TGraphErrors* gr);
  bool ProfilePRF_X(const TH2F* PRF_h, TGraphErrors* gr, TH1F* errors);

  /// Initialise PRF with expected params
  TF1* InitializePRF(const TString name);

  /// Get mean and FWHM for the histo
  Double_t GetFWHM(const TH1F* h, Double_t& mean);

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
  /// Name of the file from previous iteration
  TString _Prev_iter_name;
  /// Previous iteration output to extract PRF
  TFile*  _Prev_iter_file;

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

  /// Cluster vars
  /// Position of the cluster
  Float_t _clust_pos[Nclusters];
  /// X position of the cluster
  Float_t _x[Nclusters];
  /// Position of the "clean" cluster
  /** e.g. average 2 neighbour diagonals **/
  Float_t _cluster_av[Nclusters];
  /// X position of the evareged cluster
  Float_t _x_av[Nclusters];
  /// Position of the track
  Float_t _track_pos[Nclusters];
  /// Residuals (X_track-X_cluster)
  Float_t _residual[Nclusters];
  /// charge in the cluster
  Int_t   _charge[Nclusters];
  /// multiplicity of the cluster
  Int_t   _multiplicity[Nclusters];

  /// Pad vars
  /// X_track - X_pad --> X axis of the PRF
  Float_t _dx[Nclusters][10];
  /// Fraction of charge Q_pad / Q_cluster --> Y axis of PRF
  Float_t _qfrac[Nclusters][10];
  /// time of the pad
  Int_t   _time[Nclusters][10];
    
    // variables for the dEdx analysis
    Float_t _dEdx;
    TH1F* _hdEdx;
    Int_t _pad_charge[Nclusters][10];
    Int_t _pad_time[Nclusters][10];

    Int_t _pad_x[Nclusters][10];
    Int_t _pad_y[Nclusters][10];

    Int_t _wf_width[Nclusters][10];
    Int_t _wf_fwhm[Nclusters][10];

  /** Histograms **/
  /** Pad response function block **/
  /// PRF function from the previous step. Used for Chi2 fit
  TF1*  _PRF_function;
  /// PRF histoes
  TH2F* _PRF_histo;
  // PRF profiling graphs
  TGraphErrors* _PRF_graph;

  // TF1*    _PRF_function_2pad;
  // TF1*    _PRF_function_3pad;
  // TF1*    _PRF_function_4pad;

  TH2F* _PRF_histo_2pad;
  TH2F* _PRF_histo_3pad;
  TH2F* _PRF_histo_4pad;

  TGraphErrors* _PRF_graph_2pad;
  TGraphErrors* _PRF_graph_3pad;
  TGraphErrors* _PRF_graph_4pad;

  /// Pad response function in time
  TH2F* _PRF_time;

  /// analytical PRF time function
  TF1* _PRF_time_func;

  /// uncertainties of the time profile
  TGraphErrors* _PRF_time_error;

  TH1F* _PRF_time_e;

  /// Fitter class for the track and cluster fitting
  TrackFitCern* _fitter;

  /** Switchers **/
  /// Whether to use arc function for track fitting
  bool _do_linear_fit;
  bool _do_para_fit;
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

  /// whether to select diagonal; clusters
  bool _diagonal;

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

  TH2F* _PRF_histo_col[Nclusters];

  /// separate pad fit study
  TH1F* _Fit_quality_plots[3][Nclusters];
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
  TH1F* _x_pads = new TH1F("padX", "", 4, -0.03, 0.01);

  TH1F* _resol_col_x_scan_lim_mult[Nclusters][x_scan_bin];

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
