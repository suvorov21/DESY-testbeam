#include "SelectionBase.hxx"

#include <iostream>

bool SelectionBase::Initialize() {
  std::cout << "WARNING. The default selection is initialised. The result is always true" << std::endl;
  return true;
}

bool SelectionBase::SelectEvent(Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event& event) {
  (void)padAmpl;
  (void)event;

  return true;
}