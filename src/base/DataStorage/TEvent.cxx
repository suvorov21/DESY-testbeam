#include "TEvent.hxx"

 TEvent::TEvent(const TRawEvent& event) : TRawEvent(){
   for (const auto& hit : event.GetHits()) {
       fHitsPtrs[hit->GetCard()].emplace_back(std::make_shared<THit>(hit));
   }
 }

THitPtrVec TEvent::GetHitsInModule(const short module) {
    if (fHitsPtrs.find(module) == fHitsPtrs.end()) {
        // empty vector
        return {};
    }

    return fHitsPtrs[module];
}

void TEvent::AddHitPtr(const THitPtr &inhit) {
    fHitsPtrs[inhit->GetCard()].emplace_back(inhit);
}

void TEvent::AddPattern(const THitPtrVec &pattern) {
    if (pattern.empty())
        return;

    fPatterns[pattern.front()->GetCard()].emplace_back(pattern);
}
std::unordered_map<short, TPatternVec> TEvent::GetAllPatterns() {
    return fPatterns;
}

TPatternVec TEvent::GetPattern(const short module) {
    if (fPatterns.find(module) == fPatterns.end()) {
        // empty vector
        return {};
    }

    return fPatterns[module];
}
std::unordered_map<short, THitPtrVec> TEvent::GetAllHits() {
    return fHitsPtrs;
}

void TEvent::AddTrack(const TPatternVec &track) {
    fTracks.emplace_back(track);
}
