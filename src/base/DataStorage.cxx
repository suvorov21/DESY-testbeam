#include <DataStorage.hxx>


//// THIT



//// TTRACK
void TTrack::ResizeCols(){
  fc.resize(geom::nPadx);
}
void TTrack::ResizeRows(){
  fr.resize(geom::nPady);
}

void TTrack::AddColHit(THit* hit){
  fc[hit->GetCol()].push_back(hit);
};
void TTrack::AddRowHit(THit* hit){
  fr[hit->GetRow()].push_back(hit);
};


//// TEVENT
