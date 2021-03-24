#include "TRawEvent.hxx"

//// TRawEvent
TRawEvent::~TRawEvent() {
  for (auto hit:fUnusedHits) {
    if (hit)
      delete hit;
    hit = NULL;
  }
  for (auto hit:fUsedHits) {
    if (hit)
      delete hit;
    hit = NULL;
  }
}