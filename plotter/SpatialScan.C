#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../src/utils/Geom.hxx"
#include "../src/utils/SetT2KStyle.hxx"

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
  TCanvas c3("c3", "", 0, 630, 800, 630);
  TCanvas c4("c4", "", 800, 630, 800, 630);
  TString volt      = "360";
  TString field     = "275";
  TString peack     = "412";
  TString mag       = "02T";
  TString drift     = "530";

  TString input_prefix = "/eos/user/s/ssuvorov/DESY_testbeam/phi_merge/";
  vector<pair<TString, Int_t> > file_name_scan;

  auto scan_id = 5;
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
    case 5:
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_00" + "_iter" + TString::Itoa(Niter, 10) + ".root", 0));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_10" + "_iter" + TString::Itoa(Niter, 10) + ".root", 10));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_20" + "_iter" + TString::Itoa(Niter, 10) + ".root", 20));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_30" + "_iter" + TString::Itoa(Niter, 10) + ".root", 30));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_40" + "_iter" + TString::Itoa(Niter, 10) + ".root", 40));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_50" + "_iter" + TString::Itoa(Niter, 10) + ".root", 50));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_60" + "_iter" + TString::Itoa(Niter, 10) + ".root", 60));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_70" + "_iter" + TString::Itoa(Niter, 10) + ".root", 70));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_80" + "_iter" + TString::Itoa(Niter, 10) + ".root", 80));
      file_name_scan.push_back(make_pair(input_prefix + "phi_" + peack + "_90" + "_iter" + TString::Itoa(Niter, 10) + ".root", 90));
      scan_axis = "#phi [deg]";
      file_name += "_phi_" + peack;
      break;
  }
  file_name += ".root";
  auto out_file = new TFile(file_name, "RECREATE");

  int c_arr[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  int m_arr[9] = {20, 20, 20,204, 20, 20, 20, 20, 20};
  int l_arr[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

  TGraphErrors* resol_vs_dist     = new TGraphErrors();
  resol_vs_dist->SetName("resol");
  TGraphErrors* resol_vs_dist_e   = new TGraphErrors();
  std::vector<TGraphErrors*> resol_final;
  std::vector<TGraphErrors*> mean_final;

  TGraphErrors *resol_gr, *mean_gr;

  for (auto pair:file_name_scan) {
    TFile* f = new TFile(pair.first.Data(), "READ");
    if (!f)
      continue;

    resol_gr = (TGraphErrors*)f->Get("resol_final");
    mean_gr  = (TGraphErrors*)f->Get("mean");

    if (!resol_gr || !mean_gr)
      continue;

    resol_gr = (TGraphErrors*)resol_gr->Clone(Form("r_%i", pair.second));
    mean_gr  = (TGraphErrors*)mean_gr->Clone(Form("m_%i", pair.second));

    resol_final.push_back(resol_gr);
    mean_final.push_back(mean_gr);

    // scale to um
    for (int dot = 0; dot < (*(resol_final.end() - 1))->GetN(); dot++) {
      (*(resol_final.end() - 1))->GetY()[dot] *= 1e6;
      (*(resol_final.end() - 1))->GetEY()[dot] *= 1e6;
    }

    for (int dot = 0; dot < (*(mean_final.end() - 1))->GetN(); dot++) {
      (*(mean_final.end() - 1))->GetY()[dot] *= 1e6;
      (*(mean_final.end() - 1))->GetEY()[dot] *= 1e6;
    }

    float mean, RMS, mean_e;
    mean = GetAverage(*(resol_final.end() - 1), RMS, mean_e);
    cout << pair.second << "\t\t" << mean << "\t" << RMS << "\t" << mean_e << endl;

    resol_vs_dist->SetPoint(resol_vs_dist->GetN(), pair.second, mean);
    resol_vs_dist->SetPointError(resol_vs_dist->GetN()-1, 0, RMS);

    resol_vs_dist_e->SetPoint(resol_vs_dist_e->GetN(), pair.second, mean);
    resol_vs_dist_e->SetPointError(resol_vs_dist_e->GetN()-1, 0, mean_e);
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
  if (scan_id == 5)
    resol_vs_dist->GetYaxis()->SetRangeUser(0., 1900.);
  //resol_vs_dist->GetXaxis()->SetRangeUser(400., 600.);
  resol_vs_dist->GetXaxis()->SetTitle(scan_axis);
  resol_vs_dist->GetYaxis()->SetTitle("Resolution [#mum]");

  c2.cd();
  c2.SetGridx(1);
  c2.SetGridy(1);
  resol_final[0]->GetYaxis()->SetTitle("Resolution [#mum]");
  resol_final[0]->GetXaxis()->SetTitle("Column");
  resol_final[0]->GetYaxis()->SetRangeUser(100., 300.);
  resol_final[0]->SetLineColor(c_arr[0]);
  resol_final[0]->SetLineStyle(l_arr[0]);
  resol_final[0]->SetMarkerStyle(m_arr[0]);
  resol_final[0]->SetMarkerColor(c_arr[0]);
  resol_final[0]->Draw("ap");

  for (uint it = 1; it < resol_final.size(); ++it) {
    resol_final[it]->SetLineColor(c_arr[it]);
    resol_final[it]->SetLineStyle(l_arr[it]);
    resol_final[it]->SetMarkerStyle(m_arr[it]);
    resol_final[it]->SetMarkerColor(c_arr[it]);
    resol_final[it]->Draw("p same");
  }
  c2.BuildLegend();

  c3.cd();
  c3.SetGridx(1);
  c3.SetGridy(1);
  mean_final[0]->GetYaxis()->SetTitle("Bias [#mum]");
  mean_final[0]->GetXaxis()->SetTitle("Column");
  mean_final[0]->GetYaxis()->SetRangeUser(-300., 300.);
  mean_final[0]->SetLineColor(c_arr[0]);
  mean_final[0]->SetLineStyle(l_arr[0]);
  mean_final[0]->SetMarkerStyle(m_arr[0]);
  mean_final[0]->SetMarkerColor(c_arr[0]);
  mean_final[0]->Draw("ap");

  for (uint it = 1; it < resol_final.size(); ++it) {
    mean_final[it]->SetLineColor(c_arr[it]);
    mean_final[it]->SetLineStyle(l_arr[it]);
    mean_final[it]->SetMarkerStyle(m_arr[it]);
    mean_final[it]->SetMarkerColor(c_arr[it]);
    mean_final[it]->Draw("p same");
  }
  c3.BuildLegend();


  gPad->Modified();
  gPad->Update();
  c1.WaitPrimitive();


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