#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

/** @cond */
#include <unordered_map>
#include <array>
/** @cond */

#include "TRawEvent.hxx"
#include "THit.hxx"

using TPatternVec = std::vector<THitPtrVec>;
using broken_pads_t = std::vector<std::array<int, 3>>;

//! Class that contains the output from the reconstruction

//! Hits are divided into used and unused in the track
//! This class is used internally in the package and should NOT be used for I/O
//! Pedestals ARE subtracted in the waveforms
//! Types vocabulary:
//!   Pattern = vector of THit
//!             pattern is an association of THit based on the pattern recognition
//!   Patterns = map of vector of patterns.
//!              map of vector of vector of THit
//!              Index stands for Module. Patterns[moduleId][pattern][hit]
//!              All the patterns found in the event divided into modules
//!   Track = vector of Pattern
//!           vector of vector of vector of THit.
//!           Track[track][module][hit]
//!           Store here all the tracks after the matching of the patterns over modules
class TEvent : public TRawEvent {
 public:
    //ctor
    explicit TEvent() : TRawEvent() {}
    explicit TEvent(Int_t var) : TRawEvent(var) {}
    explicit TEvent(const TRawEvent& event);

    // getters
    std::unordered_map<short, THitPtrVec> GetAllHits();
    THitPtrVec GetHitsInModule(short module);

    TPatternVec GetPattern(short module);
    std::unordered_map<short, TPatternVec> GetAllPatterns();

    std::vector<TPatternVec> GetTracks() {return fTracks;}
    // setters

    void AddHitPtr(const THitPtr& inhit);
    void AddPattern(const THitPtrVec& pattern);

    void AddTrack(const TPatternVec& track);

// ClassDef(TEvent, 1);

 private:
    /// All the hits in event divided into modules
    std::unordered_map<short, THitPtrVec> fHitsPtrs{};

    /// Recognised patterns per module
    std::unordered_map<short, TPatternVec> fPatterns{};

    /// Patterns matched into tracks, divided into modules
    /// So the scheme is: fTracks[patternId][module][hitId]
    /// Such a structure allows independent analysis per module
    std::vector<TPatternVec> fTracks{};
};

#endif