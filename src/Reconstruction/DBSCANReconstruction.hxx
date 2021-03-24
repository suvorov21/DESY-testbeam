
#ifndef SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
#define SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_

#include "ReconstructionBase.hxx"
//! Reconstruction for passing through tracks

//! Explain it
//!
//!
//!
//!

const int MIN_DIST  = 2;
const int MIN_NODES = 1;
const int MAX_NODES = 200;
const int MAX_NODES_TOT = 300;

struct Node{
    THit* hit;
    int  c = -999;  // cluster ID
    int id = -999;  // node ID
};

struct DB_Cluster{
    int size = 0;
    int id   = -999;      // cluster ID
    // TODO make it a vector to prevent overflow
    int nodes[2000] = {-999}; // node ID of nodes belonging to cluster
};

class DBSCANReconstruction: public ReconstructionBase {
 public:
  DBSCANReconstruction();
  virtual ~DBSCANReconstruction() {;}

  virtual bool Initialize(int verbose);
  /// Main function of the reconstruction
  virtual bool SelectEvent(TEvent* event);
  /// Fill THits with maximum amplitude and time. Create Nodes
  std::vector<Node> FillNodes(TEvent* event);
  /// Merge nodes Nodes clusters
  std::vector<Node> FindClusters(std::vector<Node> raw_nodes);


  // virtual std::vector<int> FillWFs(const TEvent* event, Node n);
  virtual double MeasureDistance(Node a, Node b);

  /// Search for large enough clusters
  virtual std::vector<DB_Cluster> FindClustersLargerThan(std::vector<Node> nodes,
                                                         int minNodes
                                                         );
  /// Assotiate nodes with clusters
  virtual std::vector <Node> UpdateNodes(std::vector <DB_Cluster> clusters,
                                         std::vector <Node> nodes);
  virtual void DrawNodes(std::vector<Node> nodes);
  /// Store the output in TEven format
  virtual bool FillOutput(TEvent* event,
                          std::vector<Node> nodes,
                          std::vector<DB_Cluster> clusters);

 private:

};

#endif  // SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
