#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../utils/SetT2KStyle.hxx"

using namespace std;

const Int_t Niter = 5;

float GetAverage(TH1F* h, float& RMS);

void Spatial() {
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
  TCanvas c2("c2", "", 800, 0, 800, 630);
  TCanvas c3("c3", "", 0, 630, 800, 630);
  TCanvas c4("c4", "", 800, 630, 800, 630);

  TFile* file_in[Niter];
  TString prefix_in = "/eos/user/s/ssuvorov/DESY_testbeam/";
  TString file_name = "test";

  TGraphErrors* resol_vs_iter = new TGraphErrors();
  TGraphErrors* trackQ_vs_iter = new TGraphErrors();
  TGraphErrors* prfQ_vs_iter = new TGraphErrors();

  for (auto i = 0; i < Niter; ++i) {
    TString temp_name = prefix_in + file_name + "_iter" + TString::Itoa(i, 10) + ".root";
    file_in[i] = new TFile(temp_name, "READ");

    // fill resol vs iter
    TH1F* resol_final = (TH1F*)file_in[i]->Get("resol_final");
    resol_final->SetName(Form("resol_final_iter_%i", i));
    Float_t mean, RMS;
    mean = GetAverage(resol_final, RMS);
    resol_vs_iter->SetPoint(resol_vs_iter->GetN(), i, 1e6*mean);
    resol_vs_iter->SetPointError(resol_vs_iter->GetN()-1, 0., 1e6*RMS);

    // fill track quality vs iter
    TH1F* trackQ = (TH1F*)file_in[i]->Get("Chi2_Track");
    trackQ->SetName(Form("trackQ_iter_%i", i));
    trackQ_vs_iter->SetPoint(trackQ_vs_iter->GetN(), i, trackQ->GetMean());
    trackQ_vs_iter->SetPointError(trackQ_vs_iter->GetN()-1, 0, trackQ->GetRMS());

    // fill PRF quality vs iter
    TGraphErrors* PRF_gr = (TGraphErrors*)file_in[i]->Get("PRF_graph");
    PRF_gr->SetName(Form("PRF_graph_iter_%i", i));
    TF1* fit = PRF_gr->GetFunction("PRF_function");
    prfQ_vs_iter->SetPoint(prfQ_vs_iter->GetN(), i, fit->GetChisquare() / fit->GetNDF());
    //prfQ_vs_iter->SetPointError(resol_vs_iter->GetN()-1, i, trackQ->GetRMS());
  } // loop over iterations

  c1.cd();
  c1.SetGridx(1);
  c1.SetGridy(1);
  resol_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  resol_vs_iter->Draw("ap");

  c2.cd();
  c2.SetGridx(1);
  c2.SetGridy(1);
  trackQ_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  trackQ_vs_iter->Draw("ap");

  c3.cd();
  c3.SetGridx(1);
  c3.SetGridy(1);
  prfQ_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  prfQ_vs_iter->Draw("ap");

  c3.WaitPrimitive();
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