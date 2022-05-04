
#ifndef SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
#define SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_

#include "ReconstructionBase.hxx"

#include "TF1.h"
//! Reconstruction for passing through tracks

//! Explain it
//!
//!
//!
//!

/// distance between nodes to be associated
const int MIN_DIST = 2;
const int MIN_NODES = 1;
/// minimum number of hits for a pattern to be saved
const int MIN_HITS_PER_PATT = 15;

struct Node {
    THitPtr hit;
    int c = -999;  // cluster ID
    int id = -999;  // node ID
};

class DBSCANReconstruction : public ReconstructionBase {
 public:
    DBSCANReconstruction();

    bool Initialize(int verbose) override;
    /// Main function of the reconstruction
    bool ReconstructEvent(const std::shared_ptr<TEvent> &event) override;
    /// Fill THits with maximum amplitude and time. Create Nodes
    std::vector<Node> FillNodes(const THitPtrVec &module);
    /// Merge nodes Nodes clusters
    std::vector<Node> FindClusters(std::vector<Node> &raw_nodes);

    /// Measure distance between two nodes
    virtual double MeasureDistance(const Node &a, const Node &b);

    /// Match patterns in different modules
    void MatchModules(const std::shared_ptr<TEvent> &event);

    /// if trajectories can be fit together
    bool fitTogether(const TPattern&, const TPattern&);

    /// get adjacent modules
    std::vector<short> getNeighboursMM(ushort i);

 private:
    short chargeThresholdFit_{400};
    std::unique_ptr<TF1> trackFitFunction_;
    int fitQthreshold_{400};

};

#endif  // SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
