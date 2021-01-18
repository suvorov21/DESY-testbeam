#ifndef SRC_BASE_TRACKFITTER_HXX
#define SRC_BASE_TRACKFITTER_HXX

/** @cond */
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
/** @endcond */

#include "TEvent.hxx"
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
 * and ILC-like approach (depricated at the moment). User can use any of them
 * or specify the own approach.
 */
class TrackFitterBase {
public:
  /// shape of the track
  /// 1. Arc (circle)
  /// 2. parabolic (pol2)
  /// 3. linear (pol1)
  enum TrackShape {
    arc = 0,
    parabola,
    linear
  };

  TrackFitterBase(TrackShape shape,
                  bool invert,
                  int verbose,
                  int it);
  virtual ~TrackFitterBase() {;}

  /// set the shape of the track
  void SetTrackShape(TrackShape shape) {_shape = shape;};

  /// General function for fitting the cluster
  Double_t FitCluster();
  /// General function for fitting the whole track
  TF1* FitTrack();

  TrackFitterBase(const TrackFitterBase& fit) {
    (void)fit;
    std::cerr << "Copy constructor is depricated" << std::endl; exit(1);
  }


protected:
  // /// Fitting function for track going up
  // TF1*    _circle_function_up;
  // /// Fitting function for track going down
  // TF1*    _circle_function_dn;

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
               TH1F* _PRF_time_error,
               Float_t angle
               );
  virtual ~TrackFitCern() {;}

  /// Cluster fitter
  double FitCluster(const std::vector<THit*>& col,
                    const int& cluster,
                    const double& pos
                    );

  /// Track fitter
  TF1* FitTrack(const std::vector<TCluster*>& clusters,
                const int& miss_id = -1
                );

  /// Set array of PRFs for copmlicated patterns
  void SetPRFarr(TF1* f[], int n);

  void SetComplicatedPatternPRF(bool v) {_complicated_pattern_PRF = v;}
  void SetIndividualPRF(bool v) {_individual_column_PRF = v;}

  TrackFitCern(const TrackFitCern& fit):TrackFitterBase(fit) {(void)fit;
    std::cerr << "Copy constructor is depricated" << std::endl; exit(1);}
  bool operator==(const TrackFitCern* fit){(void)fit;
    std::cerr << "Comparison is depricated" << std::endl; exit(1);}

protected:
  /// Pad Response Function analytical description
  TF1* _PRF_function;

  /// PRF uncertainty graph
  TGraphErrors* _PRF_graph;

  /// Pad Responce Function in time
  TF1* _PRF_time_func;
  /// histogram with errors for time PRF
  TH1F* _PRF_time_error;

  /// bounds of the PRF that are reliable. Outside values are not used
  float _fit_bound;

  /// Weather to put into account uncertainty on charge with sqrt(Q)
  bool _charge_uncertainty;

  /// angle of the cluster
  Float_t _angle;

  /// Array of PRFs for complicated patterns
  TF1** _PRF_arr;

  /// PRF arr size
  int _PRF_size;

  /// Wheather to use individual PRFs for columns
  bool _individual_column_PRF;

  /// Wheather to use individual PRFs for complicated patterns e.g. 2by1 3 by1
  bool _complicated_pattern_PRF;

  /// Axis to convert track position into correction bin
  TAxis* _ax;

  /// array of corrections should be put here
  float _corr[36][50] = {{}};
};

#endif // SRC_BASE_TRACKFITTER_HXX