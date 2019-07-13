#include "SelectionBase.hxx"

bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event &event) {
  (void)padAmpl;
  (void)event;

  return true;
}