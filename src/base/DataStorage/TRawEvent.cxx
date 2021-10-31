#include "TRawEvent.hxx"

//// TRawEvent
TRawEvent::~TRawEvent() {
  for (auto hit:fHits) {
    hit.reset();
  }
  fHits.clear();
}

TRawEvent::TRawEvent(const TRawEvent* event) :
  fHits(event->fHits),
  ID(event->GetID()) {;}
