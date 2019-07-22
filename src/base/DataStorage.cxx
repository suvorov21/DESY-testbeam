#include <DataStorage.hxx>


//// THIT



//// TTRACK
TTrack::~TTrack() {
  for (auto hit:fhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }

  for (uint i = 0; i < fc.size(); ++i)
    for (auto hit:fc[i])
      hit = NULL;

  for (uint i = 0; i < fr.size(); ++i)
    for (auto hit:fr[i])
      hit = NULL;

}


//// TEVENT
TEvent::~TEvent() {
  for (auto hit:funusedhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }
  for (auto track:ftracks) {
    if (track)
      delete track;
    track = NULL;
  }
}
