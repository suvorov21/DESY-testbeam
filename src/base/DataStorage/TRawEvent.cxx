#include "TRawEvent.hxx"

TRawEvent::TRawEvent(const TRawEvent* event) :
  fHits(event->fHits),
  ID(event->GetID()) {}
