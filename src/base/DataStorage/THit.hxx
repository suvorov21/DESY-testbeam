#ifndef SRC_CLASS_THIT_HXX_
#define SRC_CLASS_THIT_HXX_

/** @cond */
#include "TROOT.h"
#include <iostream>
/** @endcond */

//! Class for storing information about each reconstructed hit.

//! Store row, column, charge, time and the whole waveform
class THit : public TObject{
 public:
  // ctor
  THit(const int col, const int row,
       const int time       = -1,
       const int q          = -1,
       const UInt_t fec     = -1,
       const UInt_t asic    = -1,
       const UInt_t channel = -1 ):
       fRow(row), fColumn(col), fTime(time), fCharge(q), fFEC(fec), fASIC(asic), fChannel(channel)  {
    fwhm      = -999;
    fw       = -999;
    for (int & i : fwf)
        fwf[i] = 0;
  }

  // default ctor
  THit() {
    fRow  = -999;
    fColumn  = -999;
    fTime  = -999;
    fCharge  = 0;

    fFEC = -1;
    fASIC = -1;
    fChannel = -1;

    for (int & i : fwf)
        fwf[i] = 0;

    fwhm = -999;
    fw = -999;
  }
  // dtor
  ~THit() override {;}
  // setters
  void SetRow(int row)    {fRow = row;}
  void SetCol(int col)    {fColumn = col;}
  void SetTime(int time)  {fTime = time;}
  void SetQ(int Q)        {fCharge = Q;}
  void SetADC(int i, int adc) {
    if (i < 0 || i > 511) {
      std::cout << "ADC index out of range!\t" << i << std::endl;
      return;
    }
    fwf[i] = adc;
  }

  void SetWidth(int val) {fw = val;}
  void SetFWHM(int val) {fwhm = val;}

  // getters
  int GetRow(bool invert = false)    const;
  int GetCol(bool invert = false)    const;
  int GetTime()               const   {return fTime;}
  int GetQ()                  const   {return fCharge;}
  int GetADC(int i) const {
    if (i >= 0 && i < 512)
      return fwf[i];
    std::cout << "ADC index out of range!\t" << i << std::endl;
    return 0;
  }
  UInt_t GetFEC()             const   {return fFEC;}
  UInt_t GetASIC()            const   {return fASIC;}
  UInt_t GetChannel()         const   {return fChannel;}

  int GetWidth() const {return fw;}
  int GetFWHM() const {return fwhm;}

  ClassDef (THit,1);

 private:
  int  fRow;
  int  fColumn;
  int  fTime;
  int  fCharge;
  int  fwf[511];
  UInt_t fFEC;
  UInt_t fASIC;
  UInt_t fChannel;

  int fwhm;
  int fw;

};

#endif
