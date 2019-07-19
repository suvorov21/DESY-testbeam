#include "ReconstructionBase.hxx"

bool ReconstructionBase::Initialize() {
  std::cout << "WARNING. The default Reconstruction is initialised. The result is always true" << std::endl;
  return true;
}

bool ReconstructionBase::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                TEvent* event) {
  (void)padAmpl;
  (void)event;

  return true;
}

std::vector<std::vector<Int_t> > ReconstructionBase::GetEmptyEvent() {
  std::vector<std::vector<Int_t> > vec;
  vec.resize(geom::nPadx);
  for (uint it = 0; it < vec.size(); ++it)
    vec[it].resize(geom::nPady, 0);

  return vec;
}
