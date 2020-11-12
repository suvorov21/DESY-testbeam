#ifndef SRC_CLASS_THIT_HXX_
#define SRC_CLASS_THIT_HXX_

/** @cond */
#include "TROOT.h"
/** @endcond */

/// Class for storing information about each reconstructed hit. Store row, column, charge and time
class THit : public TObject{
 public:
  void SetRow(int row)    {fr = row;}
  void SetCol(int col)    {fc = col;}
  void SetTime(int time)  {ft = time;}
  void SetQ(int Q)        {fq = Q;}
  void SetWidth(int w)    {fw = w;}
  void SetWHM(int whm)    {fwhm = whm;}

  int GetRow(bool invert = false)    const;
  int GetCol(bool invert = false)    const;
  int GetTime()   const   {return ft;}
  int GetQ()      const   {return fq;}
  int GetWidth()  const   {return fw;}
  int GetFWHM()   const   {return fwhm;}

  THit(int col, int row, int time, int q, int w = 0, int whm = 0){
    fr = row;
    fc = col;
    ft = time;
    fq = q;
    fw = w;
    fwhm = whm;
  }

  THit(){
    fr  = -999;
    fc  = -999;
    ft  = -999;
    fq  = 0;
  }
  virtual ~THit() {;}

  ClassDef (THit,1);

 private:
  int  fr;
  int  fc;
  int  ft;
  int  fq;
  int  fw;
  int  fwhm;
};

#endif