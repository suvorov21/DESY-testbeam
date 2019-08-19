#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "DataStorage.hxx"

namespace sel
{
  int                  GetColsMaxSep (const TTrack*);
  int                  GetColsMaxGap (const TTrack*);
  std::vector <double> GetNonZeroRows(const TTrack*);
  std::vector <double> GetNonZeroCols(const TTrack*);
  std::vector <double> GetFitParams  (const TTrack*);
  std::vector <double> Get3DFitParams(const TTrack*);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
