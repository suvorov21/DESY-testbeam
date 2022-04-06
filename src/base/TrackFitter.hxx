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

    TrackFitterBase();

    /// set the shape of the track
    void SetTrackShape(TrackShape shape) { _shape = shape; };

    void SetInversion(const bool var) { _invert = var; };
    void SetVerbosity(const int &var) { _verbose = var; };

    /// General function for fitting the cluster
    virtual std::pair<Double_t, Double_t> FitCluster(const THitPtrVec &col,
                                                     const double &pos) = 0;
    /// General function for fitting the whole track
    [[nodiscard]] virtual std::shared_ptr<TF1> FitTrack(const TClusterPtrVec &clusters,
                                                        const int &miss_id) = 0;

 protected:
    /// Arc function used for track fitting
    TF1 *_circle_function{nullptr};

    /// Fitter verbosity
    Int_t _verbose{0};

    /// GEOMETRY
    /// Either inverted geometry (rows/columns) should be used
    bool _invert{false};

    TrackShape _shape{TrackShape::parabola};
};

class TrackFitCern : public TrackFitterBase {
 public:
    /// Cluster fitter
    std::pair<Double_t, Double_t> FitCluster(const THitPtrVec &col,
                                             const double &pos
    ) override;

    /// Track fitter
    std::shared_ptr<TF1> FitTrack(const TClusterPtrVec &clusters,
                                  const int &miss_id
    ) override;

    /// Set array of PRFs for complicated patterns
    void SetPRFarr(TF1 *f[], int n);

    void SetComplicatedPatternPRF(bool v) { _complicated_pattern_PRF = v; }
    void SetIndividualPRF(bool v) { _individual_column_PRF = v; }
    void SetPRF(const TF1 *var) { _prf_function = var; }
    void SetPRFErrors(const TGraphErrors *var) { _prf_graph = var; }
    void SetPRFtimeFunc(const TF1 *var) { _prf_time_func = var; }
    void SetPRFtimeGError(TH1F *var) { _prf_time_error = var; }
    void SetFitBound(const Double_t &var) { _fit_bound = var; }
    void SetChargeUncertainty(const bool var) { _charge_uncertainty = var; }
    void SetAngle(const Double_t &var) { _angle = var; }

 protected:
    /// Pad Response Function analytical description
    const TF1 *_prf_function{nullptr};

    /// PRF uncertainty graph
    const TGraphErrors *_prf_graph{nullptr};

    /// Pad Response Function in time
    const TF1 *_prf_time_func{nullptr};
    /// histogram with errors for time PRF
    TH1F *_prf_time_error{nullptr};

    /// bounds of the PRF that are reliable. Outside values are not used
    Double_t _fit_bound{0};

    /// Weather to put into account uncertainty on charge with sqrt(Q)
    bool _charge_uncertainty{false};

    /// angle of the cluster
    Double_t _angle{0};

    /// Array of PRFs for complicated patterns
    TF1 **_prf_arr{nullptr};

    /// PRF arr size
    int _prf_size{0};

    /// Whether to use individual PRFs for columns
    bool _individual_column_PRF{false};

    /// Whether to use individual PRFs for complicated patterns e.g. 2by1 3 by1
    bool _complicated_pattern_PRF{false};
};

#endif // SRC_BASE_TRACKFITTER_HXX