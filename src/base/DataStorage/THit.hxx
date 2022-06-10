#ifndef SRC_CLASS_THIT_HXX_
#define SRC_CLASS_THIT_HXX_

/** @cond */
#include "TROOT.h"
#include <iostream>
#include <array>
/** @endcond */

#include "TRawHit.hxx"

class THit;
/// shared pointer to THit
using THitPtr = std::shared_ptr<THit>;
/// vector of shared pointers to THit, e.g. event or cluster
using THitPtrVec = std::vector<THitPtr>;

/// broken pads are arrays of [module, col, row]
using broken_pads_t = std::vector<std::array<int, 3>>;

//! Class for storing information about each reconstructed hit.

//! Store row, column, charge, time and the whole waveform
//! Pedestals ARE subtracted in the waveforms
class THit : public TRawHit {
 public:
    // ctor
    THit(const short col, const short row,
         const short card = 0, const short chip = 0, const short channel = 0,
         const short time = -1,
         const short q = -1) : TRawHit(card, chip, channel),
        fRow(row), fCol(col), fTimeMax(time), fQMax(q) {
        fwhm = -999;
        fw = -999;
    }

    // default ctor
    THit() {
        fRow = -999;
        fCol = -999;
        fTimeMax = -999;
        fQMax = 0;

        fwhm = -999;
        fw = -999;
    }

    explicit THit(const TRawHit* rhs);

    /// redefine charge and time of the hit with the maximum in given time range
    void FindMaxInTime(const int &low, const int &high);

    // setters
    void SetRow(short row) { fRow = row; }
    void SetCol(short col) { fCol = col; }
    void SetTimeMax(short time) { fTimeMax = time; }
    void SetQMax(short Q) { fQMax = Q; }

    void SetWidth(short val) { fw = val; }
    void SetFWHM(short val) { fwhm = val; }

    // getters
    short GetRow(bool invert = false) const;
    short GetCol(bool invert = false) const;
    short GetTimeMax() const { return fTimeMax; }
    short GetQMax() const { return fQMax; }

    short GetWidth() const { return fw; }
    short GetFWHM() const { return fwhm; }

// ClassDef (THit, 1);

 private:
    /// position within MM
    short fRow;
    short fCol;

    /// Time and charge of the maximum
    short fQMax;
    short fTimeMax;

    /// WF metrics
    short fwhm;
    short fw;
};

#endif
