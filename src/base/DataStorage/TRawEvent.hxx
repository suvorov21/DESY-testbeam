#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

#include "THit.hxx"

class TRawEvent : public TObject{
 public:
  //ctor
  explicit TRawEvent(){;}
  explicit TRawEvent(Int_t var): ID(var) {;}
  // dtor
  virtual ~TRawEvent();
  // getters
  Int_t GetID() const  {return ID;}
  std::vector <THit*>   GetHits()       const  {return fUsedHits;}
  std::vector <THit*>   GetUnusedHits() const  {return fUnusedHits;}
  // setters
  void SetHits(const std::vector <THit*>& inhits )        {fUsedHits = inhits;}
  void SetUnusedHits(const std::vector <THit*>& inhits)   {fUnusedHits = inhits;}
  void SetID(Int_t var) {ID = var;}

  void AddHit(THit* hit) {fUsedHits.push_back(hit);}
  void AddUnusedHit(THit* hit) {fUnusedHits.push_back(hit);}

  ClassDef (TRawEvent,1);

 private:
  /// tracks coming out of reconstruction
  std::vector <THit*>   fUsedHits;
  /// unused hits.
  std::vector <THit*>   fUnusedHits;

  Int_t ID;
};

#endif