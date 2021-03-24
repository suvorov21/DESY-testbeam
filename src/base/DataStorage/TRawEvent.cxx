#include "TRawEvent.hxx"

//// TRawEvent
TRawEvent::~TRawEvent() {
  for (auto hit:fHits) {
    if (hit)
      delete hit;
    hit = NULL;
  }
  fHits.clear();
}