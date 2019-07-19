#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
#include "CrossingReconstruction.hxx"
#include "DBSCANReconstruction.hxx"
#include "Selection.hxx"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"

/// Spatial resolution analysis
class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna(int argc, char** argv);
  virtual ~SpatialResolAna() {;}

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const Event event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();
 private:
  /// Previous iteration output to extract PRF
  TFile*  _Prev_iter_file;
  /// PRF function from the previous step. Used for Chi2 fit
  TF1*    _PRF_function;

  /// Chi2 function of the track fit
  TH1F* _Chi2_track;

  /// Residuals X_track - X_fit histoes
  TH1F* _resol_col_hist[geom::nPadx];
  TH1F* _resol_col_hist_except[geom::nPadx];

  TH1F* _resol_col_hist_2pad[geom::nPadx];
  TH1F* _resol_col_hist_2pad_except[geom::nPadx];

  TH1F* _resol_col_hist_3pad[geom::nPadx];
  TH1F* _resol_col_hist_3pad_except[geom::nPadx];

  TH1F* _residual_mean;
  TH1F* _residual_sigma;
  TH1F* _residual_sigma_2pad;
  TH1F* _residual_sigma_3pad;

  TH1F* _residual_sigma_unbiased;
  TH1F* _residual_sigma_biased;

  /// PRF histoes
  TH2F* _PRF_histo;
  TH2F* _PRF_histo_2pad;
  TH2F* _PRF_histo_3pad;
  TH2F* _PRF_histo_4pad;

  /// Average uncertainty from the previous iteration
  Float_t _uncertainty;

  // PRF profiling graphs
  TGraphErrors* _PRF_graph;

  /// vector of events IDs that passed the Reco and selection
  std::vector<Int_t> _passed_events;

  /// Chi2 track scan delta
  const float   scan_delta    = 0.1;
  /// Chi2 track scan steps
  const int     scan_Nsteps   = 1000;
  /// Chi2 track scan step
  const double  scan_step     = 2. * scan_delta / scan_Nsteps;

  // [units are meters]
  const float prf_min     = -0.018;
  const float prf_max     = 0.018;
  const int   prf_bin     = 120;

  const float resol_min   = -0.008;
  const float resol_max   = 0.008;
  const int   resol_bin   = 75.;

  const float fit_bound_left  = -0.015;
  const float fit_bound_right =  0.015;

  const float default_error   = 0.001;
  const float one_pad_error   = 0.002;
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
