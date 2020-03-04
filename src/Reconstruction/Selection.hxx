#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

/** @cond */
#include <TCanvas.h>
/** @endcond */

#include "TEvent.hxx"
#include "Geom.hxx"

namespace sel
{
  bool                 CrossingTrackSelection(const TTrack*,
                                              bool invert   = false,
                                              int _verbose  = 0);

  int                  GetColsMaxSep (const TTrack*, bool invert = false);
  int                  GetColsMaxGap (const TTrack*, bool invert = false);
  std::vector <double> GetNonZeroRows(const TTrack*, bool invert = false);
  std::vector <double> GetNonZeroCols(const TTrack*, bool invert = false);
  std::vector <double> GetFitParams  (const TTrack*, bool invert = false);
  std::vector <double> GetFitParamsXZ(const TTrack*, bool invert = false);
  std::vector <double> Get3DFitParams(const TTrack*, bool invert = false);

  // 40 ns/bin   7 cm /us  -->   0.007 cm/ns ---> 0.28 cm / bin
  static const float v_drift_est    = 0.28;
  static const float horizontal_cut = 0.1;
  static const float vertical_cut   = 0.2;
}

#endif  // SRC_SELECTION_SELECTION_HXX_
