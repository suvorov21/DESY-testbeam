#include "TEvent.hxx"

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