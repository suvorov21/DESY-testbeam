#ifndef SRC_BASE_RECONSTRUCTIONBASE_HXX_
#define SRC_BASE_RECONSTRUCTIONBASE_HXX_

/** @cond */
#include <vector>
#include <iostream>
#include <string>
/** @endcond */

#include "TROOT.h"

#include "Geom.hxx"
#include "DataStorage.hxx"


struct Hit{
  int r = -1; // row
  int c = -1; // col
  int q = 0; // charge
  int t = -1; // time
  int id = -1; // id
};

struct Track{
  std::vector<Hit> hits;             // all hits.
  std::vector<std::vector<int>>  c; // id of hits in each column
  std::vector<std::vector<int>>  r; // id of hits in each row

  void ResizeCols(){
    c.resize(geom::nPadx);
  }
  void ResizeRows(){
    r.resize(geom::nPady);
  }
  void ResizeHits(int size){
    hits.resize(size);
  }
  void PushBackCol(int col, int id){
    c[col].push_back(id);
  }
  void PushBackRow(int row, int id){
    r[row].push_back(id);
  }
  void SetHit(int id, Hit hit){
    hits[id] = hit;
  }

};

// Obsolete
/*
struct TwoD {
  Int_t A[geom::nPadx][geom::nPady] = {0};
  Int_t T[geom::nPadx][geom::nPady] = {0};
};
*/

/// The Reconstruction output structure
/** by default it's a vector of 2D event displays */
struct Event {
  std::vector <Hit> hits;       // all hits;
  std::vector <Track> tracks;  // all tracks;  ... right now strict criteria to ensure track quality.
  void ResizeTracks(int size){
    tracks.resize(size);
  }
  void ResizeHits(int size){
    hits.resize(size);
  }
  void SetHit(int id, Hit hit){
    hits[id] = hit;
  }
};

// It is important to include Selection.hxx after the stucts have been defined.
#include "Selection.hxx"

/// Template for the Reconstruction class
class ReconstructionBase {
 public:
  ReconstructionBase() {;}
  virtual ~ReconstructionBase() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                  TEvent* event);

  virtual std::vector<std::vector<Int_t>> GetEmptyEvent();
};

#endif  // SRC_BASE_RECONSTRUCTIONBASE_HXX_
