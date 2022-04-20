#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "TCluster.hxx"
#include "Geom.hxx"

const int CHARGE_THR_FOR_GAP = 0;
namespace TrackSel
{
  /// A default selection of a "nice" crossing track
  /** A cut is put on maximum multiplicity
  * and on the existence of the gaps in the cluster.
  * "Gap" is defined as a missed column or a row in the cluster, while hits
  * around this "missed pad" are detected.
  */
  bool CrossingTrackSelection(const TClusterPtrVec &track,
                            const int  &max_mult  = 10,
                            const float&max_mean_mult  = 6.,
                            const bool &cut_gap   = true,
                            const float &max_phi  = -1.,
                            const float &max_theta= -1.,
                            const std::vector<std::pair<int, int>>& _broken_pads = std::vector<std::pair<int, int>>{},
                            const bool &invert    = false,
                            const int  &verbose   = 0);

  /// Return a maximum cluster multiplicity of the track
  void GetMultiplicity (const TClusterPtrVec &track,
                       Float_t& m_mean,
                       Int_t& m_max
  );
  /// Look if there a gap in the clusters (missed row/column)
  bool GetNoGap (const TClusterPtrVec &track,
                const std::vector<std::pair<int, int>>& _broken_pads = std::vector<std::pair<int, int>>{},
                const bool &invert = false
                );

  bool GetNoGapVector(const std::vector<int>& v);

  /// Fit the track with a linear approximation
  std::vector <double> GetFitParams(const TClusterPtrVec &track,
                                    bool invert = false);

  /// Fit the track in the plane perpendicular to MM
  std::vector <double> GetFitParamsXZ(const TClusterPtrVec &track,
                                      bool invert = false);

  /// Get the track angle within the MM plane
  double GetLinearPhi(const TClusterPtrVec &track,
                      bool invert = false);

  /// Get the track angle w.r.t. MM
  double GetLinearTheta(const TClusterPtrVec &track,
                        bool invert = false);

  // 25 MHz --> 40 ns/bin   7 cm /us  -->   0.007 cm/ns ---> 0.28 cm / bin
  // 50 Mhz --> ... --> 0.14 cm / bin
  static const float v_drift_est    = 0.28;
}

#endif  // SRC_SELECTION_SELECTION_HXX_
