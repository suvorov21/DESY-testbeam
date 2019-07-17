#include "SelectionBase.hxx"

bool SelectionBase::Initialize() {
  std::cout << "WARNING. The default selection is initialised. The result is always true" << std::endl;
  return true;
}

bool SelectionBase::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                Event& event) {
  (void)padAmpl;
  (void)event;

  // 3D vector
  std::vector<std::vector<int> > twoD_temp = GetEmptyEvent();
  twoD_temp[31][35] = 500;
  event.twoD.push_back(twoD_temp);

  // vector of 2D array
  TwoD EventDisplay;
  EventDisplay.A[31][35] = 500;
  event.twoD_vector.push_back(EventDisplay);

  return true;
}

std::vector<std::vector<Int_t> > SelectionBase::GetEmptyEvent() {
  std::vector<std::vector<Int_t> > vec;
  vec.resize(geom::nPadx);
  for (uint it = 0; it < vec.size(); ++it)
    vec[it].resize(geom::nPady, 0);

  return vec;
}
