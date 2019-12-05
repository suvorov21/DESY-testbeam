// #include required here
#ifndef Geom_h
#define Geom_h

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


/*  static const int nPady = 36;
    static const int nPadx = 32;
    static const float dy = 0.011; //[m]
    static const float dx = 0.010; //[m]

    static const int pedestal = 250;

    static const float y_pos[nPady] = {-0.1925, -0.1815, -0.1705, -0.1595,
        -0.1485, -0.1375, -0.1265, -0.1155, -0.1045, -0.0935, -0.0825, -0.0715,
        -0.0605, -0.0495, -0.0385, -0.0275, -0.0165, -0.0055, 0.0055, 0.0165,
        0.0275, 0.0385, 0.0495, 0.0605, 0.0715, 0.0825, 0.0935, 0.1045, 0.1155,
        0.1265, 0.1375, 0.1485, 0.1595, 0.1705, 0.1815, 0.1925};

    static const float x_pos[nPadx] = {-0.155, -0.145, -0.135, -0.125, -0.115,
        -0.105, -0.095, -0.085, -0.075, -0.065, -0.055, -0.045, -0.035, -0.025,
        -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075,
        0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155};

    static const float MM_dy = 0.1925;
    static const float MM_dx = 0.155;
*/




    static const int Nsamples = 511;

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
