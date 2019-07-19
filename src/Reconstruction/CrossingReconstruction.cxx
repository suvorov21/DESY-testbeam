#include "CrossingReconstruction.hxx"

CrossingReconstruction::CrossingReconstruction(): ReconstructionBase() {
  ;
}

bool CrossingReconstruction::Initialize() {
  std::cout << "Initialize crossing Reconstruction............";

  std::cout << "done" << std::endl;
  return true;
}

bool CrossingReconstruction::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                TEvent* event) {
  // TODO a lot of work
  (void)padAmpl;
  (void)event;

  return true;
}
