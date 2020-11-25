#include "TCluster.hxx"

//// TCLUSTER
TCluster::~TCluster() {
  for (auto hit:fhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }

  fhits.clear();
}