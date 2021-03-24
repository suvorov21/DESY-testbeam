#include "ReconstructionBase.hxx"

bool ReconstructionBase::Initialize(int verbose) {
  (void)verbose;
  std::cout << "WARNING. The default Reconstruction is initialised. The result is always true" << std::endl;
  return true;
}

bool ReconstructionBase::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                TRawEvent* event) {
  (void)padAmpl;
  (void)event;

  return true;
}
