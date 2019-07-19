#include <DataStorage.hxx>


//// THIT



//// TTRACK
TTrack::~TTrack() {
  for (auto hit:fhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }

  for (uint i = 0; i < fc.size(); ++i)
    for (auto hit:fc[i])
      hit = NULL;

  for (uint i = 0; i < fr.size(); ++i)
    for (auto hit:fr[i])
      hit = NULL;

}

void TTrack::ResizeCols(){
  fc.resize(geom::nPadx);
}
void TTrack::ResizeRows(){
  fr.resize(geom::nPady);
}

void TTrack::AddColHit(THit* hit){
  fc[hit->GetCol()].push_back(hit);
}
void TTrack::AddRowHit(THit* hit){
  fr[hit->GetRow()].push_back(hit);
}


//// TEVENT
TEvent::~TEvent() {
  for (auto track:ftracks) {
    if (track)
      delete track;
    track = NULL;
  }
}
