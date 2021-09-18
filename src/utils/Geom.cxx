#include "Geom.hxx"

using namespace geom;

Double_t geom::GetYpos(int it_y, bool invert) {
  if ((!invert && it_y >= geom::nPady) ||
      (invert && it_y >= geom::nPadx) || it_y < 0) {
    std::cerr << "ERROR. geom::GetYpos(). Wrong Index " <<  it_y << "\t" << invert << std::endl;
    exit(1);
  }
  if (!invert)
    return geom::y_pos[it_y];
  else
    return geom::x_pos[it_y];
}

Double_t geom::GetXposPad(const THit* h, bool invert, Float_t angle) {
  return geom::GetXpos(h->GetCol(invert), invert) * TMath::Cos(angle) -
         geom::GetYpos(h->GetRow(invert), invert) * TMath::Sin(angle);
}

Double_t geom::GetYposPad(const THit* h, bool invert, Float_t angle) {
  return   geom::GetXpos(h->GetCol(invert), invert) * TMath::Sin(angle) +
           geom::GetYpos(h->GetRow(invert), invert) * TMath::Cos(angle);
}

Double_t geom::GetXpos(int it_x, bool invert) {
  if ((!invert && it_x >= geom::nPadx) ||
      (invert && it_x >= geom::nPady) || it_x < 0) {
    std::cerr << "ERROR. geom::GetXpos(). Wrong Index " <<  it_x << "\t" << invert << std::endl;
    exit(1);
  }
  if (!invert)
    return geom::x_pos[it_x];
  else
    return geom::y_pos[it_x];
}

int geom::GetNColumn(bool invert) {
  if (!invert)
    return nPadx;
  else
    return nPady;
}

int geom::GetNRow(bool invert) {
  return GetNColumn(!invert);
}
