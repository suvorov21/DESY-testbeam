#ifndef SRC_CLASS_TTRACK_HXX_
#define SRC_CLASS_TTRACK_HXX_

#include "THit.hxx"

//! Class for storing reconstructed tracks.

//! It contains vector of hits associated with this track. Also there are vectors of pointers to THit in the hit columns and rows
class TTrack : public TObject{
 public:
  void AddHit(THit* hit);

  // WARNING memory like with GetHits() observed
  std::vector<THit*> GetHits()                const     {return fhits;}
  std::vector<std::vector<THit*>> GetCols ()  const     {return fc;}
  std::vector<std::vector<THit*>> GetRows ()  const     {return fr;}

  TTrack(){;}
  virtual ~TTrack();

  ClassDef (TTrack,1);

 private:
  std::vector<THit*> fhits;               // all hits.
  std::vector<std::vector<THit*>> fc;     // id of hits in each column
  std::vector<std::vector<THit*>> fr;     // id of hits in each row

};

#endif