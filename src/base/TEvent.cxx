#include "TEvent.hxx"

TEvent::TEvent(const TRawEvent* event) {
  fUnusedHits = event->GetHits();
}