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

struct Node{
    int  x = -999;  // row
    int  y = -999;  // column
    int  t = -999;  // time 
    int  q = -999;  // charge
    int  c = -999;  // cluster ID
    int id = -999;  // node ID
};

struct Cluster{
    int size = 0;
    int id   = -999;      // cluster ID 
    int nodes[1000] = {-999}; // node ID of nodes belonging to cluster
};

class DBSCANReconstruction: public ReconstructionBase {
 public:
  DBSCANReconstruction();
  virtual ~DBSCANReconstruction() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event &event);
  virtual std::vector<Node> FillNodes(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples]);
  virtual double MeasureDistance(Node a, Node b);
  virtual std::vector<Node> FindClusters(std::vector<Node> nodes);
  virtual std::vector<Cluster> FindClustersLargerThan(std::vector<Node> nodes, int minNodes);
  virtual std::vector <Node> UpdateNodes(std::vector <Cluster> clusters, std::vector <Node> nodes);
  virtual bool CheckQuality(std::vector<Node> nodes);
  virtual void DrawNodes(std::vector<Node> nodes);
  virtual bool FillOutput(std::vector<Node> nodes, int numTracks, Event& event);

 private:

};

#endif  // SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
