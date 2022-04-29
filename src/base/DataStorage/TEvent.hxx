#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

#include "TRawEvent.hxx"
#include "THit.hxx"

//! Class that contains the output from the reconstruction

//! Hits are divided into used and unused in the track
//! This class is used internally in the package and should NOT be used for I/O
//! Pedestals ARE subtracted in the waveforms
class TEvent : public TRawEvent {
 public:
    //ctor
    explicit TEvent() : TRawEvent() {}
    explicit TEvent(Int_t var) : TRawEvent(var) {}
    explicit TEvent(const TRawEvent& event);
//     FIXME is it needed?
//    explicit TEvent(const TRawEvent* event);

    // getters
    THitPtrVec GetAllHits() const { return fHitsPtrs; }
    THitPtrVec GetUsedHits() const { return fUsedHits; }
    THitPtrVec GetUnusedHits() const { return fUnusedHits; }
    // setters
    void SetUsedHits(const THitPtrVec &inhits) { fUsedHits = inhits; }
    void SetUnusedHits(const THitPtrVec &inhits) { fUnusedHits = inhits; }
    void AddHitPtr(const THitPtr& inhit ) {fHitsPtrs.emplace_back(inhit);}

    void AddUsedHit(const THitPtr &hit) { fUsedHits.push_back(hit); }
    void AddUnusedHit(const THitPtr &hit) { fUnusedHits.push_back(hit); }

// ClassDef(TEvent, 1);

 private:
    /// tracks coming out of reconstruction
    THitPtrVec fUsedHits;
    /// unused hits.
    THitPtrVec fUnusedHits;
    /// Hit management with smart pointers
    THitPtrVec fHitsPtrs;
};

#endif