#ifndef SRC_CLASS_TRAWEVENT_HXX_
#define SRC_CLASS_TRAWEVENT_HXX_

#include "THit.hxx"

//! Class that contains pointers to all the hits in the event

//! This class is used for the I/O with a ROOT file
class TRawEvent : public TObject{
 public:
  //ctor
  explicit TRawEvent(){;}
  explicit TRawEvent(Int_t var): ID(var) {;}
  explicit TRawEvent(const TRawEvent* event) {fHits = event->fHits;}
  // dtor
  virtual ~TRawEvent();
  // getters
  Int_t GetID() const  {return ID;}
  std::vector <THit*>   GetHits()       const  {return fHits;}
  // setters
  void SetHits(const std::vector <THit*>& inhits )        {fHits = inhits;}
  void SetID(Int_t var) {ID = var;}

  void AddHit(THit* hit) {fHits.push_back(hit);}

  ClassDef (TRawEvent,1);

 protected:
  /// vector of hits in event
  std::vector <THit*> fHits;

  /// Event Id
  UInt_t ID;
};

#endif