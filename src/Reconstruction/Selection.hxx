#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "THit.hxx"
#include "TTrack.hxx"
#include "TEvent.hxx"
#include "Geom.hxx"

namespace sel
{
  int                  GetColsMaxSep (const TTrack*, bool invert = false);
  int                  GetColsMaxGap (const TTrack*, bool invert = false);
  std::vector <double> GetNonZeroRows(const TTrack*, bool invert = false);
  std::vector <double> GetNonZeroCols(const TTrack*, bool invert = false);
  std::vector <double> GetFitParams  (const TTrack*, bool invert = false);
  std::vector <double> GetFitParamsXZ(const TTrack*);
  std::vector <double> Get3DFitParams(const TTrack*, bool invert = false);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
