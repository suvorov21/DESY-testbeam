#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "TCluster.hxx"
#include "Geom.hxx"

const int CHARGE_THR_FOR_GAP = 0;
class TrackSel {
    /// phi angle of the fit
    Double_t phi_{999};
    /// theta angle of the fit
    Double_t theta_{999};

    /// function to use for a raw track fit
    std::unique_ptr<TF1> fit_;
    /// maximum multiplicity
    Double_t max_mult_{0};
    /// Max mean multiplicity
    Double_t max_mean_mult_{0};
    /// Whether to cut on gap
    bool cut_gap_{false};
    /// maximum phi angle
    Double_t max_phi_{0};
    /// maximum theta angle
    Double_t max_theta_{0};
    /// Array of broken pads
    broken_pads_t broken_pads_;

    /// Time boundaries
    Int_t time_min_;
    Int_t time_max_;
    /// whether the track is inverted
    bool invert_{false};
    /// Verbosity level
    int verbose_{0};

 public:
    TrackSel(const Double_t max_mult,
             const Double_t max_mean_mult,
             const bool cut_gap,
             const Double_t max_phi,
             const Double_t max_theta,
             broken_pads_t &broken_pads,
             const int time_min,
             const int time_max,
             const bool invert,
             const int verbose) : invert_(invert),
                                  max_mult_(max_mult),
                                  max_mean_mult_(max_mean_mult),
                                  cut_gap_(cut_gap),
                                  max_phi_(max_phi),
                                  max_theta_(max_theta),
                                  broken_pads_(broken_pads),
                                  time_min_(time_min),
                                  time_max_(time_max),
                                  verbose_(verbose) {
        fit_ = std::make_unique<TF1>("preFit", "pol2", -1., 1.);
    }

    /// Reset selection vars
    void Reset();

    ///getters
    Double_t GetPhi() const {return phi_;}
    Double_t GetTheta() const {return theta_;}

    /// A default selection of a "nice" crossing track
    /** A cut is put on maximum multiplicity
    * and on the existence of the gaps in the cluster.
    * "Gap" is defined as a missed column or a row in the cluster, while hits
    * around this "missed pad" are detected.
    */
    bool CrossingTrackSelection(const TClusterPtrVec &track);

    /// Return a maximum cluster multiplicity of the track
    static void GetMultiplicity(const TClusterPtrVec &track,
                                Float_t &m_mean,
                                Int_t &m_max
    );
    /// return min and max time for leading pad
    std::pair<int, int> GetMinMaxTime(const TClusterPtrVec &track);
    /// Look if there a gap in the clusters (missed row/column)
    bool GetNoGap(const TClusterPtrVec &track);

    static bool GetNoGapVector(const std::vector<int> &v);

    void FitXY(const TClusterPtrVec &track);
    void FitXZ(const TClusterPtrVec &track);

    /// Drift velocity estimator
    static const float v_drift_est;
};

#endif  // SRC_SELECTION_SELECTION_HXX_
