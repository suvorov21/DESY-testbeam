// #include required here
#ifndef Geom_h
#define Geom_h

/** @cond */
#include <iostream>
#include "TMath.h"
/** @endcond */

#include "THit.hxx"

//           31 ___________
//             |           |
//             |           |
//         y   |           |
//             |           |
//             |           |          __
//           0 |___________|         |__| 1.0 cm
//             0            35       1.1 cm
//                 x
//

/// Geometry definition
/** Number of pads, dimensions and time. */
namespace geom
{
    // number of pads
    static const int nPadx = 36;
    static const int nPady = 32;
    // pad size
    // (11.18+2*0.1/2 = 11.28) mm
    static const float dx = 0.01128; //[m]
    // (10.09+2*0.1/2 = 10.19) mm
    static const float dy = 0.01019; //[m]

    static const int pedestal = 250;

    static const float x_pos[nPadx] = {-0.1974, -0.18612, -0.17484, -0.16356,
        -0.15228, -0.141, -0.12972, -0.11844, -0.10716, -0.09588, -0.0846,
        -0.07332, -0.06204, -0.05076, -0.03948, -0.0282, -0.01692, -0.00564,
        0.00564, 0.01692, 0.0282, 0.03948, 0.05076, 0.06204, 0.07332, 0.0846,
        0.09588, 0.10716, 0.11844, 0.12972, 0.141, 0.15228, 0.16356, 0.17484,
        0.18612, 0.1974};

    static const float y_pos[nPady] = {-0.157945, -0.147755, -0.137565,
        -0.127375, -0.117185, -0.106995, -0.096805, -0.086615, -0.076425,
        -0.066235, -0.056045, -0.045855, -0.035665, -0.025475, -0.015285,
        -0.005095, 0.005095, 0.015285, 0.025475, 0.035665, 0.045855, 0.056045,
        0.066235, 0.076425, 0.086615, 0.096805, 0.106995, 0.117185, 0.127375,
        0.137565, 0.147755, 0.157945};

    static const int Nsamples = 511;
    static const int Nsamples_saclay = 510;

    float GetYposPad(const THit* h, bool invert = false, Float_t angle = 0);
    float GetXposPad(const THit* h, bool invert = false, Float_t angle = 0);

    float GetYpos(int it_y, bool invert = false);

    float GetXpos(int it_x, bool invert = false);

    int GetMaxColumn(bool invert = false);
    int GetMaxRow(bool invert = false);
}

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
    const float a2_inv  =   1.0652823;
    const float a3_inv  =   1.2173058;
}

#endif
