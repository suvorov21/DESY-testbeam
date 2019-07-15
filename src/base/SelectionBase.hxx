#ifndef SELECTIONBASE_H
#define SELECTIONBASE_H

#include "TROOT.h"

#include "Geom.hxx"

// the selection output structure
// by default it's a vector of 2D event displays
struct Event {
  std::vector<Int_t> twoD[geom::nPadx][geom::nPady];
};

class SelectionBase {
public:
  SelectionBase() {;}
  virtual ~SelectionBase() {;}

  bool Initialize();
  bool SelectEvent(Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event& event);

};

#endif