#ifndef SRC_CLASS_TRAWEVENT_HXX_
#define SRC_CLASS_TRAWEVENT_HXX_

#include "THit.hxx"

//! Class that contains pointers to all the hits in the event

//! This class is used for the I/O with a ROOT file
class TRawEvent : public TObject{
 public:
  //ctor
  explicit TRawEvent(): ID(0) {}
  explicit TRawEvent(Int_t var): ID(var) {}
  explicit TRawEvent(const TRawEvent* event);
  // dtor
  ~TRawEvent() override;
  // getters
  UInt_t GetID() const  {return ID;}
  std::vector <std::shared_ptr<THit>>   GetHits()       const  {return fHits;}
  // setters
  void SetHits(const std::vector<std::shared_ptr<THit>>& inhits )        {fHits = inhits;}
  virtual void SetID(Int_t var) {ID = var;}

  void AddHit(const std::shared_ptr<THit>& hit) {fHits.push_back(hit);}

  ClassDef (TRawEvent,1);

 protected:
  /// vector of hits in event
  std::vector <std::shared_ptr<THit>> fHits;

  /// Event Id
  UInt_t ID;
};

#endif