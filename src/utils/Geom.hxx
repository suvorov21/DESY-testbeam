// #include required here
#ifndef Geom_h
#define Geom_h

/** @cond */
#include <iostream>
/** @endcond */

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

    /* obsolete!!! restricted!!!
    static const int nPadx = 36;
    static const int nPady = 32;
    static const float dx = 0.011; //[m]
    static const float dy = 0.010; //[m]

    static const int pedestal = 250;

    static const float x_pos[nPadx] = {-0.1925, -0.1815, -0.1705, -0.1595,
        -0.1485, -0.1375, -0.1265, -0.1155, -0.1045, -0.0935, -0.0825, -0.0715,
        -0.0605, -0.0495, -0.0385, -0.0275, -0.0165, -0.0055, 0.0055, 0.0165,
        0.0275, 0.0385, 0.0495, 0.0605, 0.0715, 0.0825, 0.0935, 0.1045, 0.1155,
        0.1265, 0.1375, 0.1485, 0.1595, 0.1705, 0.1815, 0.1925};

    static const float y_pos[nPady] = {-0.155, -0.145, -0.135, -0.125, -0.115,
        -0.105, -0.095, -0.085, -0.075, -0.065, -0.055, -0.045, -0.035, -0.025,
        -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075,
        0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155};

    static const float MM_dx = 0.1925;
    static const float MM_dy = 0.155;
    */


    static const int Nsamples = 511;

    float GetYpos(int it_y, bool invert = false);

    float GetXpos(int it_x, bool invert = false);

    int GetMaxColumn(bool invert = false);

    // To represent geometry on canvas
    /*
    static const int wx = 420; // width of canvas in pix
    static const int wy = 340; // height of canvas in pix
    static const int times = 2;
    static const float convx = 1.;//0.85*float(wx)/(nPadx*dx); // conversion factor
    static const float convy = convx; // conversion factor
    */
}

namespace units {
    const float clight  =   299792458.;
    const float B       =   0.2;
}

#endif
