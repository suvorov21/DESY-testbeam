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
  int  fr;
  int  fc;
  int  ft;
  int  fq;

  void SetRow(int row)    {fr = row;}
  void SetCol(int col)    {fc = col;}
  void SetTime(int time)  {ft = time;}
  void SetQ(int Q)        {fq = Q;}

  int GetRow() {return fr;}
  int GetCol() {return fc;}
  int GetTime(){return ft;}
  int GetQ()   {return fq;}

  THit(){
    fr  = -999;
    fc  = -999;
    ft  = -999;
    fq  = 0;
  }
  virtual ~THit() {;}
};

class TTrack{
  public:
  std::vector<THit*> fhits; // all hits. 
  std::vector<std::vector<THit*>> fc;    // id of hits in each column
  std::vector<std::vector<THit*>> fr;    // id of hits in each row
  void ResizeCols();
  void ResizeRows();
  void AddColHit(THit* hit);
  void AddRowHit(THit* hit);
  void SetHits(std::vector<THit*> inhits)      {fhits = inhits;}
  std::vector<THit*> GetHits()                 {return fhits;}
  std::vector<THit*> GetColHits(int col)       {return fc[col];}
  std::vector<THit*> GetRowHits(int row)       {return fr[row];}
  std::vector<std::vector<THit*>> GetCols ()   {return fc;}
  std::vector<std::vector<THit*>> GetRows ()   {return fr;}

  TTrack(){
    ResizeRows();
    ResizeCols();
  }
  virtual ~TTrack() {;}
};

class TEvent{
  public:
  std::vector <THit*>   fhits;
  std::vector <TTrack*> ftracks;

  void SetHits(std::vector <THit*> inhits )       {fhits = inhits;}
  void SetTracks(std::vector <TTrack*> intracks)  {ftracks = intracks;}
  std::vector <THit*> GetHits()                   {return fhits;}
  std::vector <TTrack*> GetTracks()               {return ftracks;}

  TEvent(){;}
  virtual ~TEvent() {/*implement removal of objects...*/;}
};

#endif // SRC_CLASS_DATASTORAGE_HXX_
