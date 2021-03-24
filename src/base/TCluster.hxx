#ifndef SRC_CLASS_TCLUSTER_HXX_
#define SRC_CLASS_TCLUSTER_HXX_

#include "THit.hxx"

//! Class for storing clusters in the reconstructed tracks.

//! It contains vector of hits associated with this cluster.
class TCluster {
 public:
  void AddHit(THit* hit) {_hits.push_back(hit);};
  void SetX(Float_t x) {_x = x;}
  void SetY(Float_t y) {_y = y;}
  void SetYE(Float_t ye) {_y_error = ye;}
  void SetPos(Float_t x, Float_t y, Float_t y_e) {
    _x = x; _y = y; _y_error = y_e;
  }

  void SetCharge(Int_t q) {_charge = q;}
  void AddCharge(Int_t q) {_charge += q;}

  /// Get vector of hits
  std::vector<THit*> GetHits() const {return _hits;}
  /// Get size of the cluster == number of hits
  Ssiz_t GetSize() const {return _hits.size();}
  /// array operator
  THit* operator[](size_t n) { return _hits[n]; }

  /// loop iterator starting point
  typename std::vector<THit*>::iterator begin() {
    return _hits.begin();
  }

  /// loop iterator end
  typename std::vector<THit*>::iterator end() {
    return _hits.end();
  }

  Float_t GetX() {return _x;}
  Float_t GetY() {return _y;}
  Float_t GetYE() {return _y_error;}

  Int_t GetCharge() {return _charge;}

  TCluster(){_x = -999; _y = -999; _y_error = -999;}
  TCluster(THit* pad) {_hits.push_back(pad);}
  virtual ~TCluster();

 private:
  std::vector<THit*> _hits;               // all hits.
  Float_t _x;
  Float_t _y;
  Float_t _y_error;
  Int_t   _charge;

};

#endif