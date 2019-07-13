#ifndef SELECTIONBASE_H
#define SELECTIONBASE_H

#include "TROOT.h"

#include "../../utils/Geom.hxx"

// the selection output structure
// by default it's a vector of 2D event displays
struct Event {
  std::vector<Int_t> twoD[geom::nPadx][geom::nPady];
};

class SelectionBase {
public:
  SelectionBase() {;}
  virtual ~SelectionBase() {;}

  bool Initialise() {;}
  bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event &event);

};

#endif