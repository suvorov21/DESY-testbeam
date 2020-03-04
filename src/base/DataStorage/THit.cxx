#include "THit.hxx"

int THit::GetRow(bool invert) const {
  if (!invert)
    return fr;
  else
    return fc;
}

int THit::GetCol(bool invert) const {
  if (!invert)
    return fc;
  else
    return fr;
}