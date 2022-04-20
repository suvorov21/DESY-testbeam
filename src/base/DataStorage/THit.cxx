#include "THit.hxx"

int THit::GetRow(bool invert) const {
    if (!invert)
        return fRow;
    return fColumn;
}

int THit::GetCol(bool invert) const {
    if (!invert)
        return fColumn;
    return fRow;
}

void THit::FindMaxInTime(const int &low, const int &high) {
    for (auto time = low; time < high; ++time) {
        if (fwf[time] > fCharge) {
            fTime = time;
            fCharge = fwf[time];
        }
    } // over time
}