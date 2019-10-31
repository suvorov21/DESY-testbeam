#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../utils/Geom.hxx"
#include "../utils/SetT2KStyle.hxx"

using namespace std;

const Int_t Niter = 20;

float GetAverage(TGraphErrors* h, float& RMS);
float GetAverage(TGraphErrors* h, float& RMS, float& mean_e);

void SpatialIteration() {
  Int_t T2KstyleIndex = 2;
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
  TCanvas c4("c4", "", 0, 630, 1600, 630);

  TCanvas c5("c5", "", 1600, 0, 800, 630);
  TCanvas c6("c6", "", 1600, 630, 800, 630);
  TCanvas c7("c7", "", 2400, 0, 800, 630);
  TCanvas c8("c8", "", 2400, 630, 800, 630);
  TCanvas c9("c9", "", 2400, 670, 800, 630);


  TFile* file_in[Niter];
  TString prefix_in = "/eos/user/s/ssuvorov/DESY_testbeam/nom_v3/";
  TString file_name = "s_412_02T";

  TGraphErrors* resol_vs_iter     = new TGraphErrors();
  TGraphErrors* resol_vs_iter_sum = new TGraphErrors();
  TGraphErrors* resol_vs_iter_e   = new TGraphErrors();
  TGraphErrors* trackQ_vs_iter    = new TGraphErrors();
  TGraphErrors* trackQ_vs_iter_e  = new TGraphErrors();
  TGraphErrors* prfQ_vs_iter      = new TGraphErrors();

  TGraphErrors* mean_vs_x         = new TGraphErrors();
  TGraphErrors* sigma_vs_x        = new TGraphErrors();

  TH2F* PRF_fst, *PRF_lst;
  TGraphErrors* Fit_fst, *Fit_lst;

  for (auto i = 0; i < Niter; ++i) {
    c6.cd();
    TString temp_name = prefix_in + file_name + "_iter" + TString::Itoa(i, 10) + ".root";
    file_in[i] = new TFile(temp_name, "READ");

    // fill resol vs iter
    TGraphErrors* mean_final  = (TGraphErrors*)file_in[i]->Get("mean");
    TGraphErrors* resol_final = (TGraphErrors*)file_in[i]->Get("resol_final");
    resol_final->SetName(Form("resol_final_iter_%i", i));
    Float_t mean, RMS, mean_e;
    mean = GetAverage(resol_final, RMS, mean_e);
    resol_vs_iter->SetPoint(resol_vs_iter->GetN(), i, 1e6*mean);
    resol_vs_iter->SetPointError(resol_vs_iter->GetN()-1, 0., 1e6*RMS);
    resol_vs_iter_e->SetPoint(resol_vs_iter_e->GetN(), i, 1e6*mean);
    resol_vs_iter_e->SetPointError(resol_vs_iter_e->GetN()-1, 0., 1e6*mean_e);

    // fill track quality vs iter
    TH1F* trackQ = (TH1F*)file_in[i]->Get("Chi2_Track");
    trackQ->SetName(Form("trackQ_iter_%i", i));
    trackQ_vs_iter->SetPoint(trackQ_vs_iter->GetN(), i, trackQ->GetMean());
    trackQ_vs_iter->SetPointError(trackQ_vs_iter->GetN()-1, 0, trackQ->GetRMS());
    trackQ_vs_iter_e->SetPoint(trackQ_vs_iter_e->GetN(), i, trackQ->GetMean());
    trackQ_vs_iter_e->SetPointError(trackQ_vs_iter_e->GetN()-1, 0, trackQ->GetMeanError());

    // fill PRF quality vs iter
    TGraphErrors* PRF_gr = (TGraphErrors*)file_in[i]->Get("PRF_graph");
    PRF_gr->SetName(Form("PRF_graph_iter_%i", i));
    TF1* fit = PRF_gr->GetFunction("PRF_function");
    prfQ_vs_iter->SetPoint(prfQ_vs_iter->GetN(), i, fit->GetChisquare() / fit->GetNDF());

    if (i == 0) {
      // PRF
      Fit_fst = PRF_gr;
      PRF_fst = (TH2F*)file_in[i]->Get("PRF_histo");

      // sigma of the residual
      c5.cd();
      c5.SetGridx(1);
      c5.SetGridy(1);
      resol_final->SetMarkerColor(kBlack);
      resol_final->SetLineColor(kBlack);
      for (int dot=0;dot<resol_final->GetN();dot++) {
        resol_final->GetY()[dot] *= 1e6;
        resol_final->GetEY()[dot] *= 1e6;
      }
      resol_final->Draw("ap");
      resol_final->GetYaxis()->SetRangeUser(0., 500);
      resol_final->GetYaxis()->SetTitle("Resolution [#mum]");
      resol_final->GetXaxis()->SetTitle("Column");

      // mean of the residual
      c7.cd();
      c7.SetGridx(1);
      c7.SetGridy(1);
      mean_final->SetMarkerColor(kBlack);
      mean_final->SetLineColor(kBlack);
      for (int dot=0;dot<mean_final->GetN();dot++) {
        mean_final->GetY()[dot] *= 1e6;
        mean_final->GetEY()[dot] *= 1e6;
      }
      mean_final->Draw("ap");
      mean_final->GetYaxis()->SetRangeUser(-300., 300);
      mean_final->GetYaxis()->SetTitle("Bias [#mum]");
      mean_final->GetXaxis()->SetTitle("Column");
      c6.cd();
    }
    if (i == Niter - 1) {
      // PRF
      Fit_lst = PRF_gr;
      PRF_lst = (TH2F*)file_in[i]->Get("PRF_histo");

      // sigma of the residual
      c5.cd();
      resol_final->SetMarkerColor(kRed);
      resol_final->SetLineColor(kRed);
      for (int dot=0;dot<resol_final->GetN();dot++) {
        resol_final->GetY()[dot] *= 1e6;
        resol_final->GetEY()[dot] *= 1e6;
      }
      resol_final->Draw("p same");

      // mean of the residual
      c7.cd();
      mean_final->SetMarkerColor(kRed);
      mean_final->SetLineColor(kRed);
      for (int dot=0;dot<mean_final->GetN();dot++) {
        mean_final->GetY()[dot] *= 1e6;
        mean_final->GetEY()[dot] *= 1e6;
      }
      mean_final->Draw("p same");

      c6.cd();

      // mean of the residual vs X
      for (auto j = 0; j < 50; ++j) {
        TH1F* residual = (TH1F*)file_in[i]->Get(Form("resol_histo_Xscan_%i", j));
        if (!residual)
          continue;
        residual->Fit("gaus", "Q");
        TF1* fit_res = residual->GetFunction("gaus");
        if (!fit_res)
          continue;
        if (fit_res->GetChisquare() / fit_res->GetNDF() > 50.)
          continue;
        auto mean    = fit_res->GetParameter(1);
        auto mean_e  = fit_res->GetParError(1);
        auto sigma   = fit_res->GetParameter(2);
        auto sigma_e = fit_res->GetParError(2);
        mean_vs_x->SetPoint(mean_vs_x->GetN(), -0.035+j*(0.015+0.035)/50, 1e6*mean);
        mean_vs_x->SetPointError(mean_vs_x->GetN()-1, 0, 1e6*mean_e);

        sigma_vs_x->SetPoint(sigma_vs_x->GetN(), -0.035+j*(0.015+0.035)/50, 1e6*sigma);
        sigma_vs_x->SetPointError(sigma_vs_x->GetN()-1, 0, 1e6*sigma_e);
      }

    }

    // Other method - sum up residuals and make the average
    TH1F* residual = (TH1F*)file_in[i]->Get(Form("resol_histo_%i", 1));
    for (auto j = 2; j < geom::nPadx-1; ++j)
      if ((TH1F*)file_in[i]->Get(Form("resol_histo_%i", j)))
        residual->Add((TH1F*)file_in[i]->Get(Form("resol_histo_%i", j)));
    residual->Fit("gaus", "Q");
    TF1* fit_res = residual->GetFunction("gaus");
    auto sigma    = fit_res->GetParameter(2);
    auto sigma_e  = fit_res->GetParError(2);
    resol_vs_iter_sum->SetPoint(resol_vs_iter_sum->GetN(), i, 1e6*sigma);
    resol_vs_iter_sum->SetPointError(resol_vs_iter_sum->GetN()-1, 0, 1e6*sigma_e);
  } // loop over iterations

  c1.cd();
  c1.SetGridx(1);
  c1.SetGridy(1);
  resol_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  resol_vs_iter->GetXaxis()->SetTitle("Iteration");
  resol_vs_iter->GetYaxis()->SetRangeUser(-0., 500);
  resol_vs_iter->GetYaxis()->SetTitle("Resolution [#mum]");
  resol_vs_iter->Draw("ap >");
  resol_vs_iter_e->Draw("p same");

  c2.cd();
  c2.SetGridx(1);
  c2.SetGridy(1);
  trackQ_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  trackQ_vs_iter->GetYaxis()->SetRangeUser(0, trackQ_vs_iter->GetYaxis()->GetXmax());
  trackQ_vs_iter->Draw("ap>");
  trackQ_vs_iter_e->Draw("same p");
  c3.cd();
  c3.SetGridx(1);
  c3.SetGridy(1);
  prfQ_vs_iter->GetXaxis()->SetRangeUser(-1., Niter + 1);
  prfQ_vs_iter->Draw("ap");

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  c4.Divide(2);
  c4.cd(1);
  PRF_fst->GetYaxis()->SetRangeUser(0., 1.);
  PRF_fst->Draw("colz");
  //Fit_fst->Draw("same p");
  c4.cd(2);
  PRF_lst->GetYaxis()->SetRangeUser(0., 1.);
  PRF_lst->Draw("colz");
  //Fit_lst->Draw("same p");

  c6.cd();
  c6.SetGridx(1);
  c6.SetGridy(1);
  resol_vs_iter_sum->GetXaxis()->SetRangeUser(-1., Niter + 1);
  resol_vs_iter_sum->GetXaxis()->SetTitle("Iteration");
  resol_vs_iter_sum->GetYaxis()->SetRangeUser(-0., 500);
  resol_vs_iter_sum->GetYaxis()->SetTitle("Resolution [#mum]");
  resol_vs_iter_sum->Draw("ap");
  resol_vs_iter_sum->GetYaxis()->SetRangeUser(0., 500);

  c8.cd();
  c8.SetGridy(1);
  c8.SetGridx(1);
  mean_vs_x->GetXaxis()->SetRangeUser(-0.03, 0.01);
  mean_vs_x->GetXaxis()->SetTitle("Position [m]");
  mean_vs_x->GetYaxis()->SetRangeUser(-300., 300);
  mean_vs_x->GetYaxis()->SetTitle("Bias [#mum]");
  mean_vs_x->Draw("ap");

  c9.cd();
  c9.SetGridy(1);
  c9.SetGridx(1);
  sigma_vs_x->GetXaxis()->SetTitle("Position [m]");
  sigma_vs_x->GetXaxis()->SetRangeUser(-0.03, 0.01);
  sigma_vs_x->GetYaxis()->SetTitle("Resolution [#mum]");
  sigma_vs_x->GetYaxis()->SetRangeUser(0., 400);
  sigma_vs_x->Draw("ap");

  c4.WaitPrimitive();

  exit(0);
}

float GetAverage(TGraphErrors* h, float& RMS) {
  float mean_e;
  return GetAverage(h, RMS, mean_e);
}

float GetAverage(TGraphErrors* h, float& RMS, float& mean_e) {
  int N = 0;
  float av = 0;
  double x, y;
  for (int i = 0; i < h->GetN(); ++i) {
    h->GetPoint(i, x, y);
    av += y / h->GetN();
  }

  RMS = 0.;
  for (int i = 0; i < h->GetN(); ++i) {
    h->GetPoint(i, x, y);
    RMS += (y - av) * (y - av);
  }
  RMS = sqrt(RMS/h->GetN());
  mean_e = RMS/sqrt(h->GetN());

  return av;
}