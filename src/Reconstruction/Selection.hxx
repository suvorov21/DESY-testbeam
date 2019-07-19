#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <ReconstructionBase.hxx>

namespace sel
{
  int GetMMHits(Event event, int trackID);
  std::vector <double> GetNonZeroRows(Event event, int trackID);
  std::vector <double> GetNonZeroCols(Event event, int trackID);
  double GetFitQuality(Event event, int trackID);
}

#endif  // SRC_SELECTION_SELECTION_HXX_
