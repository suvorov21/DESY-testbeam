#ifndef SRC_CLASS_DATASTORAGE_HXX_
#define SRC_CLASS_DATASTORAGE_HXX_

/** @cond */
#include <vector>
#include <iostream>
#include <string>
/** @endcond */

#include "TROOT.h"
#include "Geom.hxx"

class THit{
 public:
  void SetRow(int row)    {fr = row;}
  void SetCol(int col)    {fc = col;}
  void SetTime(int time)  {ft = time;}
  void SetQ(int Q)        {fq = Q;}

  int GetRow()    const   {return fr;}
  int GetCol()    const   {return fc;}
  int GetTime()   const   {return ft;}
  int GetQ()      const   {return fq;}

  THit(){
    fr  = -999;
    fc  = -999;
    ft  = -999;
    fq  = 0;
  }
  virtual ~THit() {;}

 private:
  int  fr;
  int  fc;
  int  ft;
  int  fq;
};

class TTrack{
 public:
  void ResizeCols();
  void ResizeRows();
  void AddColHit(THit* hit);
  void AddRowHit(THit* hit);
  void SetHits(std::vector<THit*> inhits)               {fhits = inhits;}
  std::vector<THit*> GetHits()                const     {return fhits;}
  std::vector<THit*> GetColHits(int col)      const     {return fc[col];}
  std::vector<THit*> GetRowHits(int row)      const     {return fr[row];}
  std::vector<std::vector<THit*>> GetCols ()  const     {return fc;}
  std::vector<std::vector<THit*>> GetRows ()  const     {return fr;}

  TTrack(){
    ResizeRows();
    ResizeCols();
  }
  virtual ~TTrack() {;}

 private:
  std::vector<THit*> fhits;               // all hits.
  std::vector<std::vector<THit*>> fc;     // id of hits in each column
  std::vector<std::vector<THit*>> fr;     // id of hits in each row

};

class TEvent{
 public:
  void SetHits(std::vector <THit*> inhits )       {fhits = inhits;}
  void SetTracks(std::vector <TTrack*> intracks)  {ftracks = intracks;}
  void SetID(Int_t var) {ID = var;}
  std::vector <THit*>   GetHits()   const         {return fhits;}
  std::vector <TTrack*> GetTracks() const         {return ftracks;}
  Int_t GetID() const  {return ID;}

  TEvent(){;}
  TEvent(Int_t var): ID(var) {;}
  virtual ~TEvent() {/*implement removal of objects...*/;}

 private:
  std::vector <THit*>   fhits;
  std::vector <TTrack*> ftracks;

  Int_t ID;
};

#endif // SRC_CLASS_DATASTORAGE_HXX_
