#include "CrossingSelection.hxx"

CrossingSelection::CrossingSelection(): SelectionBase() {
  ;
}

bool CrossingSelection::Initialize() {
  std::cout << "Initialize crossing selection............";

  std::cout << "done" << std::endl;
  return true;
}

bool CrossingSelection::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                Event& event) {
  // TODO a lot of work
  (void)padAmpl;
  (void)event;

  return true;
}