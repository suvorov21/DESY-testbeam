#include "SpatialResolAna.hxx"

class InclinedTracks: public SpatialResolAna {
 public:
  InclinedTracks(int argc, char** argv);
  virtual ~InclinedTracks() {;}

  bool ProcessEvent(const TEvent* event);

  Float_t GetInclinedX(const THit* hit);
  Float_t GetInclinedY(const THit* hit);

  /// Draw the histograms of interest
  bool Draw();

  TF1* _circle_function_up;
  TF1* _circle_function_dn;
};