// #include required here
#ifndef Geom_h
#define Geom_h

/** @cond */
#include <iostream>
#include "TMath.h"
/** @endcond */

#include "TRawEvent.hxx"

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

/// Geometry definition
/** Number of pads, dimensions and time. */
class geom {
public:
  static Double_t GetYposPad(const THitPtr& h, bool invert = false, Float_t angle = 0);
  static Double_t GetXposPad(const THitPtr& h, bool invert = false, Float_t angle = 0);

  static Double_t GetYpos(int it_y, bool invert = false);

  static Double_t GetXpos(int it_x, bool invert = false);

  // time bins
  static constexpr int Nsamples = 511;
  static constexpr int Nsamples_saclay = 510;

  // number of pads
  static constexpr int nPadx = 36;
  static constexpr int nPady = 32;

  /// return max index (!!) of a column
  static int GetNColumn(bool invert = false);
  /// return a max index (!!) of a row
  static int GetNRow(bool invert = false);

private:
    geom() = default;
    // pad size
    // (11.18+2*0.1/2 = 11.28) mm
    static constexpr float dx = 0.01128; //[m]
    // (10.09+2*0.1/2 = 10.19) mm
    static constexpr float dy = 0.01019; //[m]

    static const int pedestal = 250;

    static constexpr Double_t GetXraw(int padXid);
    static constexpr Double_t GetYraw(int padYid);
};

namespace units {
    const float clight  =   299792458.;
    const float B       =   0.2;
    /** below are analytical angles
    * real angles should be corrected with the pad sizes
    const float a45     =   0.78539816;
    const float a2      =   1.1071487;
    const float a3      =   1.2490458;
    */
    const float a45     =   0.83612313;
    const float a2      =   1.1465425;
    const float a3      =   1.2783096;
    const float a32     =   1.0287271;

    const float a2_inv  =   1.0652823;
    const float a3_inv  =   1.2173058;
    const float a32_inv =   0.93503353;

    // OBSOLETE
    // const float a2      =   1.0652823;
    // const float a3      =   1.2173058;

    // const float a2_inv  =   1.1465425;
    // const float a3_inv  =   1.2783096;
}

#endif
