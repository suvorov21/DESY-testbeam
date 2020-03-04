#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "../utils/SetT2KStyle.hxx"

using namespace std;

const Int_t Niter = 20;

float GetAverage(TH1F* h, float& RMS);
float GetAverage(TH1F* h, float& RMS, float& mean_e);

void dedx() {
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
  TCanvas c4("c4", "", 800, 630, 800, 630);


  TString volt      = "360";
  TString field     = "275";
  TString peack     = "412";
  TString mag       = "02T";
  TString drift     = "430";

  TString input_prefix = "/eos/user/s/ssuvorov/DESY_testbeam/nom_v3/";

  vector<pair<TString, Float_t> > file_name_scan;
  auto scan_id = 3;
  TString scan_axis = "";
  TString file_name = input_prefix + "dedx";

  switch (scan_id) {
    case 1:
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_410" + "_dedx.root", 410));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_430" + "_dedx.root", 430));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_450" + "_dedx.root", 450));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_470" + "_dedx.root", 470));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_490" + "_dedx.root", 490));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_510" + "_dedx.root", 510));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_530" + "_dedx.root", 530));
      file_name_scan.push_back(make_pair(input_prefix + "z_"+volt+"_"+field+"_"+peack+"_"+mag+"_550" + "_dedx.root", 550));
      scan_axis = "Z position [mm]";
      file_name += "_z_" + volt + "_" + field + "_" + peack + "_" + mag;
      break;
    case 2:
      file_name_scan.push_back(make_pair(input_prefix + "g_330"+"_"+peack + "_dedx.root", 330));
      file_name_scan.push_back(make_pair(input_prefix + "g_340"+"_"+peack + "_dedx.root", 340));
      file_name_scan.push_back(make_pair(input_prefix + "g_350"+"_"+peack + "_dedx.root", 350));
      file_name_scan.push_back(make_pair(input_prefix + "g_360"+"_"+peack + "_dedx.root", 360));
      file_name_scan.push_back(make_pair(input_prefix + "g_370"+"_"+peack + "_dedx.root", 370));
      file_name += "_g_" + peack;
      scan_axis = "MM voltage [V]";
      break;
    case 3:
      file_name_scan.push_back(make_pair(input_prefix + "p_200_" + drift + "_dedx.root", 200));
      file_name_scan.push_back(make_pair(input_prefix + "p_412_" + drift + "_dedx.root", 412));
      file_name_scan.push_back(make_pair(input_prefix + "p_505_" + drift + "_dedx.root", 505));
      file_name_scan.push_back(make_pair(input_prefix + "p_116_" + drift + "_dedx.root", 116));
      file_name_scan.push_back(make_pair(input_prefix + "p_610_" + drift + "_dedx.root", 610));
      scan_axis = "Peaking time [ns]";
      file_name += "_p_" + drift;
      break;
    case 4:
      file_name_scan.push_back(make_pair(input_prefix + "m_5_" + peack + "_dedx.root", 5));
      file_name_scan.push_back(make_pair(input_prefix + "m_1_" + peack + "_dedx.root", 1));
      file_name_scan.push_back(make_pair(input_prefix + "m_2_" + peack + "_dedx.root", 2));
      file_name_scan.push_back(make_pair(input_prefix + "m_3_" + peack + "_dedx.root", 3));
      file_name_scan.push_back(make_pair(input_prefix + "m_4_" + peack + "_dedx.root", 4));
      scan_axis = "Momentum [GeV/c]";
      file_name += "_m_" + peack;
      break;
  }
  file_name += ".root";
  auto out_file = new TFile(file_name, "RECREATE");

  auto resol_vs_dist     = new TGraphErrors();
  resol_vs_dist->SetName("resol");
  auto mean_charge       = new TGraphErrors();
  mean_charge->SetName("charge");
  std::vector<TH1F*>     mult_histo;
  float max_mult = 0;

  cout << scan_axis << "\t" << "Mean charge" << "\t" << "Resolution" << "\t" << "Resolution err"  << "\t" << "Saturation" << endl;

  TH1F* resol_final;
  for (auto pair:file_name_scan) {
    auto f = new TFile(pair.first.Data(), "READ");
    resol_final = (TH1F*)f->Get("dEdx");
    resol_final->SetName(Form("dEdx_%f", pair.second));

    auto max    = resol_final->GetMaximum();
    auto max_x  = resol_final->GetBinCenter(resol_final->GetMaximumBin());
    auto FWHM   = resol_final->GetBinCenter(resol_final->FindLastBinAbove(max/2))
                - resol_final->GetBinCenter(resol_final->FindFirstBinAbove(max/2));

    resol_final->Fit("gaus", "Q", "", max_x-3*FWHM, max_x + 3*FWHM);
    TF1* fit = (TF1*)resol_final->GetFunction("gaus");
    if (!fit)
      continue;

    auto mean     = fit->GetParameter(1);
    auto sigma    = fit->GetParameter(2);
    auto mean_e   = fit->GetParError(1);
    auto sigma_e  = fit->GetParError(2);

    auto resol    = sigma / mean;
    auto resol_e  = mean_e*sigma + sigma_e*mean;
    resol_e       /= mean*mean;

    auto fst_pad    = (TH1F*)f->Get("fst_pad");
    auto saturation = 100*fst_pad->Integral(360, 1000) / fst_pad->Integral();

    TH1F* h = (TH1F*)(TH1F*)f->Get("Mult")->Clone();
    h->SetTitle(Form("%f", pair.second));
    h->SetName(Form("%f", pair.second));
    h->Scale(1/h->Integral());
    mult_histo.push_back(h);
    if (max_mult < h->GetMaximum())
      max_mult = h->GetMaximum();

    cout << pair.second << "\t\t" << mean << "\t" << resol << "\t" << resol_e << "\t" << saturation <<  endl;

    resol_vs_dist->SetPoint(resol_vs_dist->GetN(), pair.second, 100*resol);
    resol_vs_dist->SetPointError(resol_vs_dist->GetN()-1, 0, 100*resol_e);

    mean_charge->SetPoint(mean_charge->GetN(), pair.second, mean);
    mean_charge->SetPointError(mean_charge->GetN()-1, 0, mean_e);

  } // end scan

  c1.cd();
  c1.SetGridx(1);
  c1.SetGridy(1);
  resol_vs_dist->GetXaxis()->SetTitle(scan_axis);
  resol_vs_dist->GetYaxis()->SetTitle("Resolution [%]");
  resol_vs_dist->Draw("ap");
  resol_vs_dist->GetYaxis()->SetRangeUser(8., 12.);

  c2.cd();
  c2.SetGridx(1);
  c2.SetGridy(1);
  mean_charge->GetXaxis()->SetTitle(scan_axis);
  mean_charge->GetYaxis()->SetTitle("Charge [c.u.]");
  mean_charge->Draw("ap");
  if (scan_id != 2)
    mean_charge->GetYaxis()->SetRangeUser(1200., 1320.);

  c3.cd();
  resol_final->Draw();

  c4.Clear();
  c4.cd();
  mult_histo[0]->SetMaximum(max_mult*1.1);
  mult_histo[0]->SetLineColor(1);

  mult_histo[0]->Draw("hist");
  for (uint i = 1; i < mult_histo.size(); ++i) {
    mult_histo[i]->SetLineColor(i+1);
    mult_histo[i]->Draw("same hist");
  }
  c4.BuildLegend();


  gPad->Modified();
  gPad->Update();

  out_file->cd();
  resol_vs_dist->Write();
  mean_charge->Write();
  for (uint i = 0; i < mult_histo.size(); ++i) {
    mult_histo[i]->SetName(Form("mult_%i", i));
    mult_histo[i]->Write();
  }

  out_file->Close();

  c1.WaitPrimitive();

  exit(0);
}