#ifndef SRC_CLASS_TEVENT_HXX_
#define SRC_CLASS_TEVENT_HXX_

#include "TTrack.hxx"

class TEvent : public TObject{
 public:
  void SetHits(const std::vector <THit*>& inhits )       {funusedhits = inhits;}
  void SetTracks(const std::vector <TTrack*>& intracks)  {ftracks = intracks;}
  void SetID(Int_t var) {ID = var;}
  std::vector <THit*>   GetHits()   const         {return funusedhits;}
  std::vector <TTrack*> GetTracks() const         {return ftracks;}
  Int_t GetID() const  {return ID;}

  explicit TEvent(){;}
  explicit TEvent(Int_t var): ID(var) {;}
  virtual ~TEvent();

  ClassDef (TEvent,1);

 private:
  std::vector <THit*>   funusedhits; // unused hits.
  std::vector <TTrack*> ftracks;     // tracks coming out of reconstruction

  Int_t ID;
};

#endif