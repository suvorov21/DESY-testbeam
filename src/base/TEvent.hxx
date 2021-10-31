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
  explicit TEvent(Int_t var): TRawEvent(var){;}
  explicit TEvent(const TRawEvent* event) : TRawEvent(event) {;};
  // dtor
  ~TEvent() override {;}
  // getters
  std::vector <std::shared_ptr<THit>>   GetUsedHits()   const  {return fUsedHits;}
  std::vector <std::shared_ptr<THit>>   GetUnusedHits() const  {return fUnusedHits;}
  // setters
  void SetUsedHits(const std::vector <std::shared_ptr<THit>>& inhits )    {fUsedHits = inhits;}
  void SetUnusedHits(const std::vector<std::shared_ptr<THit>>& inhits)   {fUnusedHits = inhits;}
  void SetID(Int_t var) override {ID = var;}

  void AddUsedHit(std::shared_ptr<THit> hit)   {fUsedHits.push_back(hit);}
  void AddUnusedHit(std::shared_ptr<THit> hit) {fUnusedHits.push_back(hit);}

 private:
  /// tracks coming out of reconstruction
  std::vector <std::shared_ptr<THit>>   fUsedHits;
  /// unused hits.
  std::vector <std::shared_ptr<THit>>   fUnusedHits;
};

#endif