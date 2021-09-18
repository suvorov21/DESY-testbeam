#include "THit.hxx"

int THit::GetRow(bool invert) const {
    invert? return fRow: return fColumn;
  if (!invert)
    return fRow;
  return fColumn;
}

int THit::GetCol(bool invert) const {
  if (!invert)
    return fColumn;
  return fRow;
}