#include "TCluster.hxx"

//// TCLUSTER
TCluster::~TCluster() {
  for (auto hit:_hits) {
    if (hit)
      delete hit;
    hit = NULL;
  }

  _hits.clear();
}