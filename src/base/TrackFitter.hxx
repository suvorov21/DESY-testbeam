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

// define the complicated construction of pads in the column
// TODO make it as class??
typedef std::vector<std::vector<std::pair< double, std::pair<double, double> > > >  pads_t;

static pads_t v;

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
class TrackFitter {
public:
  /// All possible fitter types
  enum FitterType {
    CERN_like = 0,
    ILC_like,
    Separate_pads
  };

  enum TrackShape {
    arc = 0,
    parabola,
    linear
  };

  TrackFitter(FitterType type, TF1* func, float fit_bound,
    float uncertainty, Int_t it, Int_t verbose, bool invert,
    bool charge_uncertainty);
  virtual ~TrackFitter() {;}

  /// set the type of the fitter
  void SetType(FitterType type) {_type = type;}

  void SetTrackShape(TrackShape shape) {_shape = shape;};

  void SetDiagonal(bool d) {_diagonal = d;}

  /// General function for fitting the cluster
  Double_t FitCluster(const std::vector<THit*>& col,
    const int cluster = 0,
    double pos = 0,
    pads_t& pos_in_pad = v,
     TH1F* uncertainty = NULL);

  /// general function for fitting the whole track
  TF1* FitTrack(const std::vector<TCluster*> clusters,
    const int* mult = NULL,
    const TTrack* track = NULL,
    const double pos = 0,
    const pads_t pos_in_pad = v,
    const int miss_id = -1);

protected:
  /// Fit the whole track with CERN method
  TF1* GetTrackFitCERN(const std::vector<TCluster*> clusters,
                       const int* mult,
                       const int miss_id = -1
                       );
  /// Fit the whole track with ILC method
  TF1* GetTrackFitILC(const TTrack* track, const double pos,
                        const int miss_id = -1);
  /// Firthe whole track with independent pads
  TF1* GetTrackFitSeparatePad( const pads_t pos_in_pad, const int miss_id = -1);


  /// Extract cluster position with CERN method
  double GetClusterPosCERN(const std::vector<THit*>& col, const int cluster,
                            const double pos);
  /// Extract cluster position with ILC method
  double GetClusterPosILC(const std::vector<THit*>& col, const double pos);
  /// Fir all the pads independently with PRF
  double GetClusterPosInPad(const std::vector<THit*>& col,
    const int cluster,
    const double pos,
    pads_t& pos_in_pad,
     TH1F* uncertainty);

private:
  FitterType _type;

  TF1* _PRF_function;

  /// Fitting function for track going up
  TF1*    _circle_function_up;
  /// Fitting function for track going down
  TF1*    _circle_function_dn;

  /// bounds of the PRF that are reliable. Outside values are not used
  float _fit_bound;
  /// default value of the uncertainty (spatial resolution for prev. step)
  float _uncertainty;

  Int_t _iteration;
  Int_t _verbose;

  bool _invert;
  bool _charge_uncertainty;
  bool _do_arc_fit;
  TrackShape _shape;
  bool _diagonal;

  const float default_error   = 0.001;
  const float one_pad_error   = 0.002;

  const float sigma_pedestal = 9;
};

#endif // SRC_BASE_TRACKFITTER_HXX