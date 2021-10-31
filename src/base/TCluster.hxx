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
  auto GetHits() const {return _hits;}
  /// Get size of the cluster == number of hits
  auto GetSize() const {return _hits.size();}
  /// array operator
  THit* operator[](size_t index);

  /// loop iterator starting point
  typename std::vector<THit*>::iterator begin();

  /// loop iterator end
  typename std::vector<THit*>::iterator end();

  Float_t GetX() const {return _x;}
  Float_t GetY() const {return _y;}
  Float_t GetYE() const {return _y_error;}

  Int_t GetCharge() const {return _charge;}

  TCluster();
  explicit TCluster(THit* pad);
  virtual ~TCluster() = default;

 private:
  std::vector<THit*> _hits;               // all hits.
  Float_t _x;
  Float_t _y;
  Float_t _y_error;
  Int_t   _charge;

};

#endif