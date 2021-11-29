#include "Geom.hxx"

Double_t geom::GetYpos(int it_y, bool invert) {
  if ((!invert && it_y >= geom::nPady) ||
      (invert && it_y >= geom::nPadx) || it_y < 0) {
    TString msg = __func__;
    msg += "Wrong Index " + std::to_string(it_y) + "\t" + std::to_string(invert);
    throw std::logic_error(msg);
  }
  if (!invert)
    return geom::GetYraw(it_y);
  else
    return geom::GetXraw(it_y);
}

Double_t geom::GetXposPad(const THitPtr& h, bool invert, Float_t angle) {
  return geom::GetXpos(h->GetCol(invert), invert) * TMath::Cos(angle) -
         geom::GetYpos(h->GetRow(invert), invert) * TMath::Sin(angle);
}

Double_t geom::GetYposPad(const THitPtr& h, bool invert, Float_t angle) {
  return   geom::GetXpos(h->GetCol(invert), invert) * TMath::Sin(angle) +
           geom::GetYpos(h->GetRow(invert), invert) * TMath::Cos(angle);
}

Double_t geom::GetXpos(int it_x, bool invert) {
  if ((!invert && it_x >= geom::nPadx) ||
      (invert && it_x >= geom::nPady) || it_x < 0) {
    TString msg = __func__;
    msg += "Wrong Index " + std::to_string(it_x) + "\t" + std::to_string(invert);
    throw std::logic_error(msg);
  }
  if (!invert)
    return geom::GetXraw(it_x);
  else
    return geom::GetYraw(it_x);
}

int geom::GetNColumn(bool invert) {
  if (!invert)
    return nPadx;
  else
    return nPady;
}

int geom::GetNRow(bool invert) {
  return geom::GetNColumn(!invert);
}

constexpr Double_t geom::GetXraw(int padXid) {
  return dx * (padXid - 1.*(nPadx-1) / 2);
}

constexpr Double_t geom::GetYraw(int padYid) {
  return dy * (padYid - 1.*(nPady-1) / 2);
}
