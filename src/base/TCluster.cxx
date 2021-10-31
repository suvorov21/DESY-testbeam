#include "TCluster.hxx"
#include "Geom.hxx"

//// TCLUSTER
TCluster::TCluster(): _x(-999), _y(-999), _y_error(-999), _charge(0) {}

TCluster::TCluster(THit* pad) {
  _hits.push_back(pad);
  _x = (float_t)geom::GetXposPad(pad);
  _y = (float_t)geom::GetYposPad(pad);
  _y_error = 0;
  _charge = 0;
}

typename std::vector<THit*>::iterator TCluster::begin() {
  return _hits.begin();
}

typename std::vector<THit*>::iterator TCluster::end() {
  return _hits.end();
}

THit* TCluster::operator[](size_t index) {
  return _hits[index];
}