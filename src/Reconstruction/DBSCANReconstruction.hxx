
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

struct Node{
    int  x = -999;  // row
    int  y = -999;  // column
    int  t = -999;  // time
    int  q = -999;  // charge
    int  c = -999;  // cluster ID
    int id = -999;  // node ID
    int  w = -999;  // WF width
    int whm = -999; // WF WHM
    //std::vector<std::pair<int, int>> wf_v; //WF

};



struct Cluster{
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
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], TEvent* event);
  virtual std::vector<Node> FillNodes(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples]);
  //virtual std::vector<std::vector<int>> FillWFs(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples]);
  virtual std::vector<int> FillWFs(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Node n);
  virtual double MeasureDistance(Node a, Node b);
  virtual std::vector<Node> FindClusters(std::vector<Node> nodes);
  virtual std::vector<Cluster> FindClustersLargerThan(std::vector<Node> nodes, int minNodes);
  virtual std::vector <Node> UpdateNodes(std::vector <Cluster> clusters, std::vector <Node> nodes);
  virtual void DrawNodes(std::vector<Node> nodes);
  virtual bool FillOutput(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], std::vector<Node> nodes, std::vector<Cluster> clusters, TEvent* event);

 private:

};

#endif  // SRC_RECONSTRUCTION_DBSCANRECONSTRUCTION_HXX_
