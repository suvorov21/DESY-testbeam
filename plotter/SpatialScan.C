#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../utils/SetT2KStyle.hxx"

using namespace std;

const Int_t Niter = 0;

float GetAverage(TH1F* h, float& RMS);

void SpatialScan() {
  Int_t T2KstyleIndex = 3;
  // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
  TString localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = T2KstyleIndex;

  TStyle* t2kstyle = T2K().SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());
  gROOT->ForceStyle();

  TCanvas c1("c1", "", 0, 0, 800, 630);


  TString input_prefix = "/eos/user/s/ssuvorov/DESY_testbeam/v6/";
  vector<pair<TString, Float_t> > file_name_scan;

  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_410_iter" + TString::Itoa(Niter, 10) +  ".root", 410));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_430_iter" + TString::Itoa(Niter, 10) +  ".root", 430));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_450_iter" + TString::Itoa(Niter, 10) +  ".root", 450));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_470_iter" + TString::Itoa(Niter, 10) +  ".root", 470));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_490_iter" + TString::Itoa(Niter, 10) +  ".root", 490));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_510_iter" + TString::Itoa(Niter, 10) +  ".root", 510));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_530_iter" + TString::Itoa(Niter, 10) +  ".root", 530));
  file_name_scan.push_back(make_pair(input_prefix + "z_360_275_200_02T_550_iter" + TString::Itoa(Niter, 10) +  ".root", 550));

  TGraphErrors* resol_vs_dist = new TGraphErrors();

  for (auto pair:file_name_scan) {
    TFile* f = new TFile(pair.first.Data(), "READ");
    TH1F* resol_final = (TH1F*)f->Get("resol_final");
    resol_final->SetName(Form("resol_final_%f", pair.second));
    float mean, RMS;
    mean = GetAverage(resol_final, RMS);

    resol_vs_dist->SetPoint(resol_vs_dist->GetN(), pair.second, 1.e6*mean);
    resol_vs_dist->SetPointError(resol_vs_dist->GetN()-1, 1.e6*RMS);
  }

  c1.cd();
  resol_vs_dist->Draw("ap");
  resol_vs_dist->GetYaxis()->SetRangeUser(0., 500.);
  resol_vs_dist->GetXaxis()->SetRangeUser(400., 600.);
  gPad->Modified();
  gPad->Update();
  c1.WaitPrimitive();

}

float GetAverage(TH1F* h, float& RMS) {
  int N = 0;
  float av = 0;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    if (!h->GetBinContent(i))
      continue;
    ++N;
    av += h->GetBinContent(i);
  }
  av /= N;

  RMS = 0.;
  for (int i = 2; i <= h->GetXaxis()->GetNbins() - 1; ++i) {
    RMS += (h->GetBinContent(i) - av) * (h->GetBinContent(i) - av);
  }
  RMS = sqrt(RMS/N);

  return av;
}