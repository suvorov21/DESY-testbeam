#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "TEvent.hxx"
#include "TCluster.hxx"
#include "Geom.hxx"

namespace sel
{
  /// A default selection of a "nice" crossing track
  /** A cut is put on maximum multiplicity
  * and on the exictance of the gaps in the cluster.
  * "Gap" is defined as a missed column or a row in the cluster, while hits
  * around this "missed pad" are detected.
  */
  bool                 CrossingTrackSelection(const std::vector<TCluster*> &track,
                                              const int  &max_mult  = 6,
                                              const bool &cut_gap   = true,
                                              const float &max_phi  = -1.,
                                              const float &max_theta= -1.,
                                              const bool &invert    = false,
                                              const int  &verbose   = 0);

  /// Return a maximum cluster multiplicity of the track
  int               GetMaxMultiplicity (const std::vector<TCluster*> &track);
  /// Look if there a gap in the clusters (missed row/column)
  bool                  GetNoGap       (const std::vector<TCluster*> &track,
                                        const bool &invert = false);

  /// Fit the track with a linear approximation
  std::vector <double> GetFitParams  (const std::vector<TCluster*> &track,
                                      bool invert = false);

  /// Fit the track in the plane perpendicular to MM
  std::vector <double> GetFitParamsXZ(const std::vector<TCluster*> &track,
                                      bool invert = false);

  /// Get the track angle within the MM plane
  float GetLinearPhi(const std::vector<TCluster*> &track,
                     bool invert = false);

  /// Get the track angle w.r.t. MM
  float GetLinearTheta(const std::vector<TCluster*> &track,
                       bool invert = false);

  // 25 MHz --> 40 ns/bin   7 cm /us  -->   0.007 cm/ns ---> 0.28 cm / bin
  // 50 Mhz --> ... --> 0.14 cm / bin
  static const float v_drift_est    = 0.28;
  static const float horizontal_cut = 0.1;
  static const float vertical_cut   = 0.1;
}

#endif  // SRC_SELECTION_SELECTION_HXX_
