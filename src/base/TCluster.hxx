#ifndef SRC_CLASS_TCLUSTER_HXX_
#define SRC_CLASS_TCLUSTER_HXX_

#include "THit.hxx"

class TCluster;
using TClusterPtr = std::shared_ptr<TCluster>;
using TClusterPtrVec = std::vector<TClusterPtr>;

//! Class for storing clusters in the reconstructed tracks.

//! It contains vector of hits associated with this cluster.
class TCluster {
 public:
    void AddHit(const THitPtr &hit) { _hits.push_back(hit); };
    void SetX(Double_t x) { _x = x; }
    void SetY(Double_t y) { _y = y; }
    void SetYE(Double_t ye) { _y_error = ye; }
    void SetPos(Double_t x, Double_t y, Double_t y_e) {
        _x = x;
        _y = y;
        _y_error = y_e;
    }

    void SetCharge(Int_t q) { _charge = q; }
    void AddCharge(Int_t q) { _charge += q; }

    /// Get vector of hits
    auto GetHits() const { return _hits; }
    /// Get size of the cluster == number of hits
    auto GetSize() const { return _hits.size(); }
    /// array operator
    THitPtr operator[](size_t index);

    /// loop iterator starting point
    typename THitPtrVec::iterator begin();

    /// loop iterator end
    typename THitPtrVec::iterator end();

    Double_t GetX() const { return _x; }
    Double_t GetY() const { return _y; }
    Double_t GetYE() const { return _y_error; }

    Int_t GetCharge() const { return _charge; }

    TCluster();
    explicit TCluster(const std::shared_ptr<THit> &pad);

 private:
    /// vector pf all the hits in the cluster
    THitPtrVec _hits{};
    Double_t _x;
    Double_t _y;
    Double_t _y_error;
    Int_t _charge;

};


enum class clusterType {
    kRowColumn,
    kDiagonal,
    k3by1,
    k2by1
};

/// Class that keeps rules for track clusterisation
class Clustering {
    /// angle of a reference frame rotation
    Double_t _angle;
    /// Number of pads in a row
    int _n_pads;
    /// whether to invert MM
    bool _invert;
    /// Slope coefficient. 0 corresponds to columns/rows. 1 to diagonals and so on
    Double_t _coeff;
 public:
    Clustering(clusterType type, bool invert);

    // Split track into clusters
    /** Extract the vector of clusters from the whole track.
    * the logic of clusterisation is given with the function of the Clustering object
    * The function takes (row, column) and return a value
    * that is constant for a given cluster.
    * For example for clustering with columns the rule column == const is constant.
    * For diagonals column - row = const and so on.
    */
    TClusterPtrVec ClusterTrack(const THitPtrVec &tr) const;

    int getNpads() const {return _n_pads;}
    Double_t getAngle() const {return _angle;}

    /// Function of row and column that is constant for a given cluster
    [[nodiscard]] int GetConstant(int row, int col) const {
        if (_n_pads == 0)
            return col;
        else {
            return floor(_coeff * col - row);
        }
    }
};

#endif