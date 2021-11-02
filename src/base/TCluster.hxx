#ifndef SRC_CLASS_TCLUSTER_HXX_
#define SRC_CLASS_TCLUSTER_HXX_

#include "TRawEvent.hxx"

class TCluster;
using TClusterPtr = std::shared_ptr<TCluster>;
using TClusterPtrVec = std::vector<TClusterPtr>;

//! Class for storing clusters in the reconstructed tracks.

//! It contains vector of hits associated with this cluster.
class TCluster {
 public:
  void AddHit(const THitPtr& hit) {_hits.push_back(hit);};
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
  THitPtr operator[](size_t index);

  /// loop iterator starting point
  typename THitPtrVec::iterator begin();

  /// loop iterator end
  typename THitPtrVec::iterator end();

  Float_t GetX() const {return _x;}
  Float_t GetY() const {return _y;}
  Float_t GetYE() const {return _y_error;}

  Int_t GetCharge() const {return _charge;}

  TCluster();
  explicit TCluster(std::shared_ptr<THit> pad);
  virtual ~TCluster() = default;

 private:
  /// vector pf all the hits in the cluster
  THitPtrVec _hits;
  Float_t _x;
  Float_t _y;
  Float_t _y_error;
  Int_t   _charge;

};

#endif