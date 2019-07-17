#include "SelectionBase.hxx"

bool SelectionBase::Initialize() {
  std::cout << "WARNING. The default selection is initialised. The result is always true" << std::endl;
  return true;
}

bool SelectionBase::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                Event& event) {
  (void)padAmpl;
  (void)event;

  std::vector<std::vector<int> > twoD_temp = GetEmptyEvent();

  event.twoD.push_back(twoD_temp);
  return true;
}

std::vector<std::vector<Int_t> > SelectionBase::GetEmptyEvent() {
  std::vector<std::vector<Int_t> > vec;
  vec.resize(geom::nPadx);
  for (uint it = 0; it < vec.size(); ++it)
    vec.resize(geom::nPady);

  return vec;
}
