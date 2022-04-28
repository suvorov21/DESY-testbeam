#include "THit.hxx"

THit::THit(const TRawHit* rhs) : TRawHit(rhs) {
    fQMax = 0;
    const auto& wf = GetADCvector();
    for (auto t = 0; t < wf.size(); ++t) {
        // subtract pedestals
        auto Q = wf[t] - Geom::pedestal;
        SetADCunit(num::cast<short>(GetTime() + t),num::cast<short>(Q));
        if (Q > fQMax) {
            fQMax = num::cast<short>(Q);
            fTimeMax = num::cast<short>(GetTime() + t);
        }
    }

    auto pos = Geom::get().GetEleToPad(GetChip(),
                                       GetChannel());
    fRow = pos.second;
    fCol = pos.first;

    // width characteristics
    fwhm = 0;
    fw = 0;
    for (short t : GetADCvector()) {
        if (t > 0)
            ++fw;
        if (t > fQMax / 2)
            ++fwhm;
    }
}

short THit::GetRow(bool invert) const {
    if (!invert)
        return fRow;
    return fCol;
}

short THit::GetCol(bool invert) const {
    if (!invert)
        return fCol;
    return fRow;
}

void THit::FindMaxInTime(const int &low, const int &high) {
    fQMax = -999;
    for (auto time = low; time < high; ++time) {
        if (GetADCvector()[time] > fQMax) {
            fTimeMax = num::cast<short>(time + GetTime());
            fQMax = GetADCvector()[time];
        }
    } // over time
}