#ifndef SRC_CLASS_TCLUSTER_HXX_
#define SRC_CLASS_TCLUSTER_HXX_

#include "THit.hxx"

//! Class for storing clusters in the reconstructed tracks.

//! It contains vector of hits associated with this cluster.
class TCluster : public TObject{
 public:
  void AddHit(THit* hit) {fhits.push_back(hit);};
  void SetX(Float_t x) {_x = x;}
  void SetY(Float_t y) {_y = y;}
  void SetYE(Float_t ye) {_y_error = ye;}
  void SetPos(Float_t x, Float_t y, Float_t y_e) {
    _x = x; _y = y; _y_error = y_e;
  }

  void SetCharge(Int_t q) {_charge = q;}

  // THit& operator[] (unsigned i) {return *fhits[i];}

  std::vector<THit*> GetHits()                const     {return fhits;}

  Float_t GetX() {return _x;}
  Float_t GetY() {return _y;}
  Float_t GetYE() {return _y_error;}

  Int_t GetCharge() {return _charge;}

  TCluster(){_x = -999; _y = -999; _y_error = -999;}
  TCluster(THit* pad) {fhits.push_back(pad);}
  virtual ~TCluster();

  ClassDef (TCluster,1);

 private:
  std::vector<THit*> fhits;               // all hits.
  Float_t _x;
  Float_t _y;
  Float_t _y_error;
  Int_t   _charge;

};

#endif