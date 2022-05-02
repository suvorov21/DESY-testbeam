//
// Created by SERGEY SUVOROV on 09/04/2022.
//

#ifndef DESY_TESTBEAM_SRC_RECONSTRUCTION_PADSELECTION_HXX_
#define DESY_TESTBEAM_SRC_RECONSTRUCTION_PADSELECTION_HXX_

#include "TEvent.hxx"
#include "TCluster.hxx"

namespace PadSelection {
/// Process a cluster and return only pads that are suggested to be robust
/** E.g. function can return only 2 pads in a column.
 * Another use case is to omit pads with wrong timestamps.
 * Any user defined selection may be applied.
 */
THitPtrVec GetRobustPadsInCluster(THitPtrVec col,
                                  const broken_pads_t& broken_pads);
/// Return only robust clusters
/** E.g. apply a truncation - omit clusters with relatively large charge
 * Or put a strong upper limit on cluster charge.
 * Any condition can be specified.
 */
TClusterPtrVec GetRobustClusters(TClusterPtrVec &tr);
};

#endif //DESY_TESTBEAM_SRC_RECONSTRUCTION_PADSELECTION_HXX_
