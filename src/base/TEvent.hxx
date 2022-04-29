#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

#include "TRawEvent.hxx"

//! Class that contains the output from the reconstruction

//! Hits are divided into used and unused in the track
//! This class is used internally in the package and should NOT be used for I/O
class TEvent : public TRawEvent {
 public:
    //ctor
    explicit TEvent() : TRawEvent() {}
    explicit TEvent(Int_t var) : TRawEvent(var) {}
    explicit TEvent(const TRawEvent *event) : TRawEvent(event) {}
    explicit TEvent(const TRawEvent &event) : TRawEvent(event) {}

    // getters
    THitPtrVec GetUsedHits() const { return fUsedHits; }
    THitPtrVec GetUnusedHits() const { return fUnusedHits; }
    // setters
    void SetUsedHits(const THitPtrVec &inhits) { fUsedHits = inhits; }
    void SetUnusedHits(const THitPtrVec &inhits) { fUnusedHits = inhits; }
    void SetID(Int_t var) override { ID = var; }

    void AddUsedHit(const THitPtr &hit) { fUsedHits.push_back(hit); }
    void AddUnusedHit(const THitPtr &hit) { fUnusedHits.push_back(hit); }

 private:
    /// tracks coming out of reconstruction
    THitPtrVec fUsedHits;
    /// unused hits.
    THitPtrVec fUnusedHits;
};

#endif