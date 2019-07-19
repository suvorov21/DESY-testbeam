#include "Selection.hxx"

int sel::GetMMHits(Event event, int trackID){
  int hits = 0;
  // for(int i=0; i<geom::nPadx;i++)
  //   for(int j=0; j<geom::nPady;j++)
  //     if(event.twoD[trackID][i][j]) hits++;
  return hits;
}

std::vector <double> sel::GetNonZeroRows(Event event, int trackID){
  std::vector <double> rows;
  int rowQ = 0;
  for(auto row:event.tracks[trackID].r) if(row.size()){
    int rowQ = 0;
    for(auto q:row) rowQ+=event.tracks[trackID].hits[q].q;
    if(rowQ) rows.push_back(rowQ);
  }
  return rows;
}

std::vector <double> sel::GetNonZeroCols(Event event, int trackID){
  std::vector <double> cols;
  int colQ = 0;
  for(auto col:event.tracks[trackID].c) if(col.size()){
    int colQ = 0;
    for(auto q:col) colQ+=event.tracks[trackID].hits[q].q;
    if(colQ) cols.push_back(colQ);
  }
  return cols;
}

double sel::GetFitQuality(Event event, int trackID){
  
  TCanvas c;
  c.cd();
  TH2F    *MM      = new TH2F("MM","MM",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  
  // for(int i=0; i<geom::nPadx;i++)
  //   for(int j=0; j<geom::nPady;j++)
  //     if(event.twoD[trackID][i][j]) MM->Fill(i,j,event.twoD[trackID][i][j]);
  for(auto h:event.tracks[trackID].hits) if(h.q) MM->Fill(h.c,h.r,h.q);

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  double quality = 1.0e10;
  if (fit){
    quality = fit->GetChisquare() / fit->GetNDF();
    double k = fit->GetParameter(1);
    double b = fit->GetParameter(0);
  }

  delete MM;
  return quality;
}  
