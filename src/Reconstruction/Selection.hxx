#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <ReconstructionBase.hxx>

namespace sel
{
  int GetMMHits(const TEvent* event, int trackID);
  std::vector <double> GetNonZeroRows(const TEvent* event, int trackID);
  std::vector <double> GetNonZeroCols(const TEvent* event, int trackID);
  double GetFitQuality(const TEvent* event, int trackID);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
