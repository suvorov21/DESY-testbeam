#include "THit.hxx"

int THit::GetRow(bool invert) const {
  if (!invert)
    return fRow;
  else
    return fColumn;
}

int THit::GetCol(bool invert) const {
  if (!invert)
    return fColumn;
  else
    return fRow;
}