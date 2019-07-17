#ifndef SRC_BASE_SELECTIONBASE_HXX_
#define SRC_BASE_SELECTIONBASE_HXX_

/** @cond */
#include <vector>
#include <iostream>
#include <string>
/** @endcond */

#include "TROOT.h"

#include "Geom.hxx"

struct TwoDdisplay {
  Int_t twoD[geom::nPadx][geom::nPady];
};

/// The selection output structure
/** by default it's a vector of 2D event displays */
struct Event {
  std::vector<std::vector<std::vector<Int_t> > > twoD;
  std::vector<TwoDdisplay> twoD_vector;
};

/// Template for the selection class
class SelectionBase {
 public:
  SelectionBase() {;}
  virtual ~SelectionBase() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                  Event& event);

  virtual std::vector<std::vector<Int_t> > GetEmptyEvent();
};

#endif  // SRC_BASE_SELECTIONBASE_HXX_
