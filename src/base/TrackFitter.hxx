#ifndef SRC_BASE_TRACKFITTER_HXX
#define SRC_BASE_TRACKFITTER_HXX

/** @cond */
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
/** @endcond */

#include "Geom.hxx"
#include "TCluster.hxx"

/**
 * @brief      This class describes a track fitter.
 * The class provides two main functions for fitting both individual clusters
 * (reconstructed dot) and whole tracks (true track position).
 * Double_t FitCluster() returns position of the dot from the cluster
 * TF1* FitTrack() returns analytical function that describes the true track.
 *
 * Different methods are given with FitterType: CERN-like, individual pad fit
 * and ILC-like approach (deprecated at the moment). User can use any of them
 * or specify the own approach.
 */
class TrackFitterBase {
public:
  /// shape of the track
  /// 1. Arc (circle)
  /// 2. parabolic (pol2)
  /// 3. linear (pol1)
  enum class TrackShape {
    arc = 0,
    parabola,
    linear
  };

  TrackFitterBase(TrackShape shape,
                  bool invert,
                  int verbose,
                  int it);

  /// set the shape of the track
  void SetTrackShape(TrackShape shape) {_shape = shape;};

  /// General function for fitting the cluster
  virtual Double_t FitCluster();
  /// General function for fitting the whole track
  virtual std::shared_ptr<TF1> FitTrack();

protected:
  /// Arc function used for track fitting
  TF1* _circle_function;

  /// iteration of the fit
  Int_t _iteration;
  /// Fitter verbosity
  Int_t _verbose;

  /// GEOMETRY
  /// Either inverted geometry (rows/columns) should be used
  bool _invert;

  TrackShape _shape;

  const float default_error   = 0.001;
  const float one_pad_error   = 0.002;

  // const float sigma_pedestal = 9;
};

class TrackFitCern: public TrackFitterBase {
public:

  TrackFitCern(TrackShape shape,
               bool invert,
               int verbose,
               int it,
               TF1* PRF,
               TGraphErrors* PRF_gr,
               float fit_bound,
               bool charge_uncertainty,
               TF1* PRF_time_func,
               TH1F* PRF_time_error,
               Float_t angle
               );

  TrackFitCern() : TrackFitCern(TrackShape::arc, false, 0, 0,
                                nullptr, nullptr, 0.0,
                                false, nullptr,
                                nullptr, 0.) {};

  /// Cluster fitter
  Double_t FitCluster(const THitPtrVec& col,
                     const double& pos = 0
                     );

  /// Track fitter
  std::shared_ptr<TF1> FitTrack(const TClusterPtrVec& clusters,
                                const int& miss_id = -1
                                );

  /// Set array of PRFs for complicated patterns
  void SetPRFarr(TF1* f[], int n);

  void SetComplicatedPatternPRF(bool v) {_complicated_pattern_PRF = v;}
  void SetIndividualPRF(bool v) {_individual_column_PRF = v;}

protected:
  /// Pad Response Function analytical description
  TF1* _prf_function{nullptr};

  /// PRF uncertainty graph
  TGraphErrors* _prf_graph{nullptr};

  /// Pad Response Function in time
  TF1* _prf_time_func{nullptr};
  /// histogram with errors for time PRF
  TH1F* _prf_time_error{nullptr};

  /// bounds of the PRF that are reliable. Outside values are not used
  float _fit_bound{0};

  /// Weather to put into account uncertainty on charge with sqrt(Q)
  bool _charge_uncertainty{false};

  /// angle of the cluster
  Float_t _angle{0};

  /// Array of PRFs for complicated patterns
  TF1** _prf_arr{nullptr};

  /// PRF arr size
  int _prf_size{0};

  /// Whether to use individual PRFs for columns
  bool _individual_column_PRF{false};

  /// Whether to use individual PRFs for complicated patterns e.g. 2by1 3 by1
  bool _complicated_pattern_PRF{false};

  /// Axis to convert track position into correction bin
  TAxis* _ax{nullptr};

  /// array of corrections should be put here
  float _corr[36][50] = {{}};
};

#endif // SRC_BASE_TRACKFITTER_HXX