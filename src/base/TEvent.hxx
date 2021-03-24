#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

#include "TRawEvent.hxx"

//! Class that contains the output from the reconstruction

//! Hits are devided into used and unused in the track
//! This class is used internally in the package and should NOT be used for I/O
class TEvent: public TRawEvent {
 public:
  //ctor
  explicit TEvent() : TRawEvent() {;}
  explicit TEvent(UInt_t var): TRawEvent(var){;}
  explicit TEvent(const TRawEvent* event) : TRawEvent(event) {;};
  // dtor
  virtual ~TEvent() {;};
  // getters
  std::vector <THit*>   GetUsedHits()   const  {return fUsedHits;}
  std::vector <THit*>   GetUnusedHits() const  {return fUnusedHits;}
  // setters
  void SetUsedHits(const std::vector <THit*>& inhits )    {fUsedHits = inhits;}
  void SetUnusedHits(const std::vector <THit*>& inhits)   {fUnusedHits = inhits;}
  void SetID(Int_t var) {ID = var;}

  void AddUsedHit(THit* hit)   {fUsedHits.push_back(hit);}
  void AddUnusedHit(THit* hit) {fUnusedHits.push_back(hit);}

 private:
  /// tracks coming out of reconstruction
  std::vector <THit*>   fUsedHits;
  /// unused hits.
  std::vector <THit*>   fUnusedHits;


};

#endif