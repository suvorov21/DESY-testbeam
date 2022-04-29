#ifndef Geom_h
#define Geom_h

/** @cond */
#include <iostream>
#include <unordered_map>

#include "TMath.h"
/** @endcond */

#include "THit.hxx"

//           31 ___________
//             |           |
//             |           |
//         y   |           |
//             |           |
//             |           |          __
//           0 |___________|         |__| 1.019 cm
//             0            35       1.128 cm
//                 x
//

class THit;
using THitPtr = std::shared_ptr<THit>;

/// Geometry definition
/** Number of pads, dimensions and time. */
class Geom {
 public:
    // Singleton stuff
    Geom(const Geom &) = delete;
    Geom &operator=(const Geom &) = delete;
    static Geom &get();
    /// extract Y position from pad
    static Double_t GetYposPad(const THitPtr &h,
                               bool invert = false,
                               Double_t angle = 0);

    /// extract X position from pad
    static Double_t GetXposPad(const THitPtr &h,
                               bool invert = false,
                               Double_t angle = 0);

    /// Get Y position of the column
    static Double_t GetYpos(int it_y, bool invert = false);

    /// Get X position of the row
    static Double_t GetXpos(int it_x, bool invert = false);

    // time bins
    static constexpr short Nsamples = 511;
    static constexpr short Nsamples_saclay = 510;

    // number of pads
    static constexpr short nPadx = 36;
    static constexpr short nPady = 32;

    static constexpr short nChip = 16;
    static constexpr short nChannel = 79;

    // pad size
    // (11.18+2*0.1/2 = 11.28) mm
    static constexpr float dx = 0.01128; //[m]
    // (10.09+2*0.1/2 = 10.19) mm
    static constexpr float dy = 0.01019; //[m]

    /// Pedestals values
    static const short pedestal = 250;

    /// return max index (!!) of a column
    static int GetNColumn(bool invert = false);
    /// return a max index (!!) of a row
    static int GetNRow(bool invert = false);

    std::pair<short, short> GetPadToEle(int x, int y);
    std::pair<short, short> GetEleToPad(int chip, int channel);

 private:
    Geom();

    /// private method to access to raw X position array
    static constexpr Double_t GetXraw(int padXid);
    /// private method to access to raw Y position array
    static constexpr Double_t GetYraw(int padYid);

    std::pair<short, short> padToEle[nPadx * nPady];
    std::pair<short, short> eleToPad[nChip * nChannel];
};

namespace units {
const float clight = 299792458.;
const float B = 0.2;

const float a45 = 0.83612313;
const float a2 = 1.1465425;
const float a3 = 1.2783096;

const float a2_inv = 1.0652823;
const float a3_inv = 1.2173058;
}

namespace num {
/// safe cast for numerical values
template<class T, class TT>
T cast(TT tt) {
    T t = static_cast<T>(tt);
    if (static_cast<TT>(t) != tt)
        throw std::bad_cast();

    if ((std::is_signed<T>::value ^ std::is_signed<TT>::value) && \
          (t < T(0) || tt < TT(0)))
        throw std::bad_cast();

    return t;
}
}

#endif
