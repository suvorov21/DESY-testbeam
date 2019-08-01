#include "TF1.h"
#include "TH2F.h"

#include "Selection.hxx"

int sel::GetColsMaxSep(const TTrack* track){
  int maxsep = 0;
  for(auto col:track->GetCols()) if(col.size()){
    int first = 0;
    int last = 0;
    std::sort(col.begin(), col.end(), [](THit* hit1, THit* hit2){return hit1->GetRow() < hit2->GetRow();});
    int diff = (*(col.end()-1))->GetRow() - (*col.begin())->GetRow();
    if (diff > maxsep) maxsep = diff;
  }
  return maxsep+1;
}

int sel::GetColsMaxGap(const TTrack* track){
  int maxgap  = 0;
  for(auto col:track->GetCols()) if(col.size()){
    std::sort(col.begin(), col.end(), [](THit* hit1, THit* hit2){return hit1->GetRow() < hit2->GetRow();});
    for (uint padID = 1; padID < col.size(); ++padID) {
      if (col[padID]->GetRow() - col[padID-1]->GetRow()-1 > maxgap)
        maxgap = col[padID]->GetRow() - col[padID-1]->GetRow()-1;

    }
  }
  return maxgap;
}

std::vector <double> sel::GetNonZeroRows(const TTrack* track){
  std::vector <double> rows;
  for(auto row:track->GetRows()) if(row.size()){
    int rowQ = 0;
    for(auto h:row) rowQ+=h->GetQ();
    if(rowQ) rows.push_back(rowQ);
  }
  return rows;
}

std::vector <double> sel::GetNonZeroCols(const TTrack* track){
  std::vector <double> cols;
  for(auto col:track->GetCols()) if(col.size()){
    int colQ = 0;
    for(auto h:col) colQ+=h->GetQ();
    if(colQ) cols.push_back(colQ);
  }
  return cols;
}

std::vector <double> sel::GetFitParams(const TTrack* track){
  std::vector <double> params;

  TH2F    *MM      = new TH2F("MM","MM",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  for(auto h:track->GetHits()) if(h->GetQ()) MM->Fill(h->GetCol(),h->GetRow(),h->GetQ());

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  double quality = 1.0e10;
  if (fit){
    quality = fit->GetChisquare() / fit->GetNDF();
    double b = fit->GetParameter(0);
    double k = fit->GetParameter(1);
    params.push_back(quality);
    params.push_back(b);
    params.push_back(k);
  }

  delete MM;
  return params;
}

//// 3D FITTING:    ---- WARNING ---- UNDER DEVELOPMENT, DO NOT USE IT!

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TH3.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TNtuple.h>
#include <Fit/Fitter.h>

using namespace ROOT::Math;

// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
  // a parametric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = t;
}

// calculate distance line-point
double distToFit(double x,double y,double z, std::vector <double> p) {
  // distance line point is D= | (xp-x0) cross  ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0)
  XYZVector xp(x,y,z);
  XYZVector x0(p[0], p[2], 0. );
  XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
  XYZVector u = (x1-x0).Unit();
  double d2 = ((xp-x0).Cross(u)).Mag2();
  return d2;
}


bool first = true;

// function Object to be minimized
struct SumDistance2 {
  const TTrack *fTrack;

  SumDistance2(const TTrack *g) : fTrack(g) {}

  // calculate distance line-point
  double distance2(double x,double y,double z, const double *p) {
     // distance line point is D= | (xp-x0) cross  ux |
     // where ux is direction of line and x0 is a point in the line (like t = 0)
     XYZVector xp(x,y,z);
     XYZVector x0(p[0], p[2], 0. );
     XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
     XYZVector u = (x1-x0).Unit();
     double d2 = ((xp-x0).Cross(u)).Mag2();
     return d2;
  }

  // implementation of the function to be minimized
  double operator() (const double *par) {
     assert(fTrack != 0);
     double sum = 0;
     for (auto h:fTrack->GetHits()) {
        double d = distance2(h->GetCol(),h->GetRow(),h->GetTime(),par);
        sum += d;
     }
     first = false;
     return sum;
  }

};

std::vector<double> sel::Get3DFitParams(const TTrack* track)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.05);

  ROOT::Fit::Fitter  fitter;

  // make the functor objet
  SumDistance2 sdist(track);
  ROOT::Math::Functor fcn(sdist,4);
  // set the function and the initial parameter values
  double pStart[4] = {1,1,1,1};
  fitter.SetFCN(fcn,pStart);
  // set step sizes different than default ones (0.3 times parameter values)
  for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

  std::vector <double> p_result;

  bool ok = fitter.FitFCN();
  if (!ok) {
     Error("line3Dfit","Line3D Fit failed");
     return p_result;
  }

  const ROOT::Fit::FitResult & result = fitter.Result();

  // std::cout << "Total final distance square: " << result.MinFcnValue() << std::endl;
  // std::cout << "Ave Dist: " << result.MinFcnValue() / (gr->GetN()) << std::endl;
  // result.Print(std::cout);

  //if(result.MinFcnValue() / (track->GetHits().size()) > 1.5) return p_result;

  // get fit parameters
  const double * parFit = result.GetParams();

  // draw the fitted line
  int n = 1000;
  double t0 = 0;
  double dt = 1000;
  TPolyLine3D *l = new TPolyLine3D(n);
  for (int i = 0; i <n;++i) {
     double t = t0+ dt*i/n;
     double x,y,z;
     line(t,parFit,x,y,z);
     l->SetPoint(i,x,y,z);
  }
  l->SetLineColor(kRed);

  double p0 = parFit[0];
  double p1 = parFit[1];
  double p2 = parFit[2];
  double p3 = parFit[3];

  p_result.push_back(p0);
  p_result.push_back(p1);
  p_result.push_back(p2);
  p_result.push_back(p3);

  bool DrawThis = true;


  std::cout << "here..." << std::endl;
  if(DrawThis){
    std::cout << "drawing..." << std::endl;
    TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");
    for (auto h:track->GetHits()){
      if(!h->GetQ()) continue;
      event3D->Fill(h->GetTime(),h->GetRow(),h->GetCol(),h->GetQ());
    }
    TCanvas *canv = new TCanvas("canv", "canv", 800, 600, 800, 600);
    canv->cd(1);
    event3D->Draw("x:y:z:c","","box2");
    TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetLimits(0,geom::nPadx);
    htemp->GetYaxis()->SetLimits(0,geom::nPady);
    htemp->GetZaxis()->SetLimits(0,500);
    htemp->SetTitle("");
    l->Draw("same");
    canv->Update();
    canv->WaitPrimitive();
    delete htemp;
    delete canv;
    delete event3D;
  }

  return p_result;
}



