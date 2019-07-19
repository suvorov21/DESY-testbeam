#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <ReconstructionBase.hxx>

namespace sel
{
  int GetMMHits(TEvent* event, int trackID);
  std::vector <double> GetNonZeroRows(TEvent* event, int trackID);
  std::vector <double> GetNonZeroCols(TEvent* event, int trackID);
  double GetFitQuality(TEvent* event, int trackID);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
