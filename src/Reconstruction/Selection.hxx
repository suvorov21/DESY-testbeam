#ifndef SRC_SELECTION_SELECTION_HXX_
#define SRC_SELECTION_SELECTION_HXX_

#include "ReconstructionBase.hxx"

namespace sel
{

  int GetMMHits(Event event, int trackID){
    int hits = 0;
    for(int i=0; i<geom::nPadx;i++)
      for(int j=0; j<geom::nPadx;j++)
        if(event.twoD[trackID][i][j]) hits++;
    return hits;
  }

  int GetNonZeroRows(Event event, int trackID){
    int rows = 0;
    for(int i=0; i<geom::nPadx;i++){
      int rowQ = 0;
      for(int j=0; j<geom::nPadx;j++){
        if(event.twoD[trackID][i][j]) rowQ+=event.twoD[trackID][i][j];
      }
      if(rowQ) rows++;
    }
    return rows;
  }

  int GetNonZeroCols(Event event, int trackID){
    int cols = 0;
    for(int i=0; i<geom::nPadx;i++){
      int colQ = 0;
      for(int j=0; j<geom::nPadx;j++){
        if(event.twoD[trackID][j][i]) colQ+=event.twoD[trackID][j][i];
      }
      if(colQ) cols++;
    }
    return cols;
  }
  
}

#endif  // SRC_SELECTION_SELECTION_HXX_
