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

const Int_t Niter = 19;

float GetAverage(TGraphErrors* h, float& RMS);
float GetAverage(TGraphErrors* h, float& RMS, float& mean_e);

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
  TCanvas c2("c2", "", 800, 0, 800, 630);
  TString volt      = "360";
  TString field     = "275";
  TString peack     = "200";
  TString mag       = "02T";
  TString drift     = "530";

  TString input_prefix = "/eos/user/s/ssuvorov/DESY_testbeam/nom_v3/";
  vector<pair<TString, Float_t> > file_name_scan;

  auto scan_id = 4;
  TString scan_axis = "";
  TString file_name = input_prefix + "SR";

  switch (scan_id) {
    case 1:
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_410_iter" + TString::Itoa(Niter, 10) +  ".root", 410));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_430_iter" + TString::Itoa(Niter, 10) +  ".root", 430));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_450_iter" + TString::Itoa(Niter, 10) +  ".root", 450));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_470_iter" + TString::Itoa(Niter, 10) +  ".root", 470));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_490_iter" + TString::Itoa(Niter, 10) +  ".root", 490));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_510_iter" + TString::Itoa(Niter, 10) +  ".root", 510));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_530_iter" + TString::Itoa(Niter, 10) +  ".root", 530));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_550_iter" + TString::Itoa(Niter, 10) +  ".root", 550));
      scan_axis = "Z position [mm]";
      file_name += "_z_" + volt + "_" + field + "_" + peack + "_" + mag;
      break;
    case 2:
      file_name_scan.push_back(make_pair(input_prefix + "g_330"+"_"+peack+"_iter" + TString::Itoa(Niter, 10) +  ".root", 330));
      file_name_scan.push_back(make_pair(input_prefix + "g_340"+"_"+peack+"_iter" + TString::Itoa(Niter, 10) +  ".root", 340));
      file_name_scan.push_back(make_pair(input_prefix + "g_350"+"_"+peack+"_iter" + TString::Itoa(Niter, 10) +  ".root", 350));
      file_name_scan.push_back(make_pair(input_prefix + "g_360"+"_"+peack+"_iter" + TString::Itoa(Niter, 10) +  ".root", 360));
      file_name_scan.push_back(make_pair(input_prefix + "g_370"+"_"+peack+"_iter" + TString::Itoa(Niter, 10) +  ".root", 370));
      scan_axis = "MM voltage [V]";
      file_name += "_g_" + peack;
      break;
    case 3:
      file_name_scan.push_back(make_pair(input_prefix + "p_200_" + drift + "_iter" + TString::Itoa(Niter, 10) +  ".root", 200));
      file_name_scan.push_back(make_pair(input_prefix + "p_412_" + drift + "_iter" + TString::Itoa(Niter, 10) +  ".root", 412));
      file_name_scan.push_back(make_pair(input_prefix + "p_505_" + drift + "_iter" + TString::Itoa(Niter, 10) +  ".root", 505));
      file_name_scan.push_back(make_pair(input_prefix + "p_116_" + drift + "_iter" + TString::Itoa(Niter, 10) +  ".root", 116));
      file_name_scan.push_back(make_pair(input_prefix + "p_610_" + drift + "_iter" + TString::Itoa(Niter, 10) +  ".root", 610));
      scan_axis = "Peaking time [ns]";
      file_name += "_p_" + drift;
      break;
    case 4:
      file_name_scan.push_back(make_pair(input_prefix + "m_5_" + peack + "_iter" + TString::Itoa(Niter, 10) + ".root", 5));
      file_name_scan.push_back(make_pair(input_prefix + "m_1_" + peack + "_iter" + TString::Itoa(Niter, 10) + ".root", 1));
      file_name_scan.push_back(make_pair(input_prefix + "m_2_" + peack + "_iter" + TString::Itoa(Niter, 10) + ".root", 2));
      file_name_scan.push_back(make_pair(input_prefix + "m_3_" + peack + "_iter" + TString::Itoa(Niter, 10) + ".root", 3));
      file_name_scan.push_back(make_pair(input_prefix + "m_4_" + peack + "_iter" + TString::Itoa(Niter, 10) + ".root", 4));
      scan_axis = "Momentum [GeV/c]";
      file_name += "_m_" + peack;
      break;
  }
  file_name += ".root";
  auto out_file = new TFile(file_name, "RECREATE");

  TGraphErrors* resol_vs_dist     = new TGraphErrors();
  resol_vs_dist->SetName("resol");
  TGraphErrors* resol_vs_dist_e   = new TGraphErrors();

  for (auto pair:file_name_scan) {
    TFile* f = new TFile(pair.first.Data(), "READ");
    auto resol_final = (TGraphErrors*)f->Get("resol_final");
    resol_final->SetName(Form("resol_final_%f", pair.second));
    float mean, RMS, mean_e;
    mean = GetAverage(resol_final, RMS, mean_e);
    cout << pair.second << "\t\t" << 1.e6*mean << "\t" << 1.e6*RMS << "\t" << 1.e6*mean_e << endl;

    resol_vs_dist->SetPoint(resol_vs_dist->GetN(), pair.second, 1.e6*mean);
    resol_vs_dist->SetPointError(resol_vs_dist->GetN()-1, 0, 1.e6*RMS);

    resol_vs_dist_e->SetPoint(resol_vs_dist_e->GetN(), pair.second, 1.e6*mean);
    resol_vs_dist_e->SetPointError(resol_vs_dist_e->GetN()-1, 0, 1.e6*mean_e);
/*

    // Other method - sum up residuals and make the average
    TH1F* residual = (TH1F*)f->Get(Form("resol_histo_%i", 1));
    for (auto i = 2; i < geom::nPadx-1; ++i)
      if ((TH1F*)f->Get(Form("resol_histo_%i", i)))
        residual->Add((TH1F*)f->Get(Form("resol_histo_%i", i)));
    residual->Fit("gaus");*/

  }


  resol_vs_dist->GetXaxis()->SetTitle(scan_axis);
  resol_vs_dist->GetYaxis()->SetTitle("Resolution [#mum]");


  c1.cd();
  c1.SetGridx(1);
  c1.SetGridy(1);
  resol_vs_dist->Draw("ap>");
  resol_vs_dist_e->Draw("p same");

  resol_vs_dist->GetYaxis()->SetRangeUser(0., 400.);
  //resol_vs_dist->GetXaxis()->SetRangeUser(400., 600.);
  resol_vs_dist->GetXaxis()->SetTitle(scan_axis);
  resol_vs_dist->GetYaxis()->SetTitle("Resolution [#mum]");
  gPad->Modified();
  gPad->Update();
  //c1.WaitPrimitive();

  out_file->cd();
  resol_vs_dist->Write();
  out_file->Close();

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