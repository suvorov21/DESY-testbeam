#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "DataStorage.hxx"

namespace sel
{
  int GetMMHits(const TEvent* event, int trackID);
  std::vector <double> GetNonZeroRows(const TEvent* event, int trackID);
  std::vector <double> GetNonZeroCols(const TEvent* event, int trackID);
  std::vector <double> GetFitParams(const TEvent* event, int trackID);
  std::vector <double> Get3DFitParams(TTrack* track);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
