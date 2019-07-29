#include <DataStorage.hxx>


//// THIT



//// TTRACK
TTrack::~TTrack() {
  for (auto hit:fhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }

  fc.clear();
  fr.clear();
}

void TTrack::AddHit(THit* hit) {
  // fill general vector
  fhits.push_back(hit);

  // fill column vector
  auto col_found = false;
  for (uint colID = 0; colID < fc.size(); ++colID) {
    if (fc[colID][0]->GetCol() == hit->GetCol()) {
      fc[colID].push_back(hit);
      col_found = true;
      break;
    }
  }
  if (!col_found) {
    std::vector<THit*> temp_v;
    temp_v.push_back(hit);
    fc.push_back(temp_v);
  }
  // fill row vector
  auto row_found = false;
  for (uint rowID = 0; rowID < fr.size(); ++rowID) {
    if (fr[rowID][0]->GetRow() == hit->GetRow()) {
      fr[rowID].push_back(hit);
      row_found = true;
      break;
    }
  }
  if (!row_found) {
    std::vector<THit*> temp_v;
    temp_v.push_back(hit);
    fr.push_back(temp_v);
  }
}


//// TEVENT
TEvent::~TEvent() {
  for (auto hit:funusedhits) {
    if (hit)
      delete hit;
    hit = NULL;
  }
  for (auto track:ftracks) {
    if (track)
      delete track;
    track = NULL;
  }
}
