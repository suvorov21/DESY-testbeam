#include "TEvent.hxx"

 TEvent::TEvent(const TRawEvent& event) : TRawEvent(){
   for (const auto& hit : event.GetHits()) {
       fHitsPtrs.emplace_back(std::make_shared<THit>(hit));
   }
 }