// #include required here
#ifndef Geom_h
#define Geom_h

//           31 ___________
//             |           |
//             |           |
//         y   |           |
//             |           |
//             |           |          __
//           0 |___________|         |__| 1.1 cm
//             0            35       1.0 cm
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

#endif
