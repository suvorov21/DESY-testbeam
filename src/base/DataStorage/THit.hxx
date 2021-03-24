#ifndef SRC_CLASS_THIT_HXX_
#define SRC_CLASS_THIT_HXX_

/** @cond */
#include "TROOT.h"
/** @endcond */

/// Class for storing information about each reconstructed hit. Store row, column, charge and time
class THit : public TObject{
 public:
  // ctor
  THit(const int col, const int row,
       const int time       = -1,
      const int q           = -1,
       const std::vector<int> wf = std::vector<int>(),
       const UInt_t fec     = -1,
       const UInt_t asic    = -1,
       const UInt_t channel = -1 ){
    fRow      =  row;
    fColumn   = col;
    fTime     = time;
    fCharge   = q;
    fwf       = wf;
    fFEC      = fec;
    fASIC     = asic;
    fChannel  = channel;
  }

  // default ctor
  THit(){
    fRow  = -999;
    fColumn  = -999;
    fTime  = -999;
    fCharge  = 0;
    fwf.clear();
  }
  // dtor
  // Assume std::vector is destroyed automatically
  virtual ~THit() {;}
  // setters
  void SetRow(int row)    {fRow = row;}
  void SetCol(int col)    {fColumn = col;}
  void SetTime(int time)  {fTime = time;}
  void SetQ(int Q)        {fCharge = Q;}
  void SetWF_v(std::vector<int> wf){fwf = wf;}

  // getters
  int GetRow(bool invert = false)    const;
  int GetCol(bool invert = false)    const;
  int GetTime()               const   {return fTime;}
  int GetQ()                  const   {return fCharge;}
  std::vector<int> GetWF_v()  const   {return fwf;}
  UInt_t GetFEC()             const   {return fFEC;}
  UInt_t GetASIC()            const   {return fASIC;}
  UInt_t GetChannel()         const   {return fChannel;}

  ClassDef (THit,1);

 private:
  int  fRow;
  int  fColumn;
  int  fTime;
  int  fCharge;
  std::vector<int> fwf;
  UInt_t fFEC;
  UInt_t fASIC;
  UInt_t fChannel;

};

#endif
