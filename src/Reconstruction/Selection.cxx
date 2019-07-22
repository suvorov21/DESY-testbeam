#include "Selection.hxx"

int sel::GetMMHits(const TEvent* event, int trackID){
  return event->GetHits().size();
}

std::vector <double> sel::GetNonZeroRows(const TEvent* event, int trackID){
  std::vector <double> rows;
  for(auto row:event->GetTracks()[trackID]->GetRows()) if(row.size()){
    int rowQ = 0;
    for(auto h:row) rowQ+=h->GetQ();
    if(rowQ) rows.push_back(rowQ);
  }
  return rows;
}

std::vector <double> sel::GetNonZeroCols(const TEvent* event, int trackID){
  std::vector <double> cols;
  for(auto col:event->GetTracks()[trackID]->GetCols()) if(col.size()){
    int colQ = 0;
    for(auto h:col) colQ+=h->GetQ();
    if(colQ) cols.push_back(colQ);
  }
  return cols;
}

double sel::GetFitQuality(const TEvent* event, int trackID){

  TH2F    *MM      = new TH2F("MM","MM",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  for(auto h:event->GetTracks()[trackID]->GetHits()) if(h->GetQ()) MM->Fill(h->GetCol(),h->GetRow(),h->GetQ());

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  double quality = 1.0e10;
  if (fit){
    quality = fit->GetChisquare() / fit->GetNDF();
    // double k = fit->GetParameter(1);
    // double b = fit->GetParameter(0);
  }

  delete MM;
  return quality;
}
