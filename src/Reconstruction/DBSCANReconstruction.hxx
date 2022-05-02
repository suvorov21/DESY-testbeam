
#ifndef SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
#define SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_

#include "ReconstructionBase.hxx"
//! Reconstruction for passing through tracks

//! Explain it
//!
//!
//!
//!

const int MIN_DIST = 2;
const int MIN_NODES = 1;
const int MIN_NODES_PER_PATT = 15;

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

    // virtual std::vector<int> FillWFs(const TEvent* event, Node n);
    virtual double MeasureDistance(const Node &a, const Node &b);

    void MatchModules(const std::shared_ptr<TEvent> &event);

 private:

};

#endif  // SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
