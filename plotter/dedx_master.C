#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TLegend.h"

#include "../utils/SetT2KStyle.hxx"

using namespace std;

const Int_t Niter = 20;

float GetAverage(TH1F* h, float& RMS);
float GetAverage(TH1F* h, float& RMS, float& mean_e);

void dedx_master() {
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
  c1.SetGridy(1);
  c1.SetGridx(1);

  TString input_prefix  = "/eos/user/s/ssuvorov/DESY_testbeam/nom_v3/";
  TString file_name     = input_prefix + "/dedx.root";
  auto out_file         = new TFile(file_name, "RECREATE");

  std::vector< std::vector< std::pair<TString, TString> > > list;

  list.push_back(std::vector<std::pair<TString, TString> > {std::make_pair("dedx_z_360_275_200_02T.root", "200 ns 0.2 T"),
                                                          std::make_pair("dedx_z_360_275_200_0T.root", "200 ns 0. T"),
                                                          std::make_pair("dedx_z_360_275_412_02T.root", "412 ns 0.2 T"),
                                                          std::make_pair("dedx_z_360_275_412_0T.root", "412 ns 0. T.")});

  list.push_back(std::vector<std::pair<TString, TString> > {std::make_pair("dedx_g_200.root", "200 ns"),
                                                          std::make_pair("dedx_g_412.root", "412 ns")});

  list.push_back(std::vector<std::pair<TString, TString> > {std::make_pair("dedx_m_200.root", "200 ns"),
                                                          std::make_pair("dedx_m_412.root", "412 ns")});

  list.push_back(std::vector<std::pair<TString, TString> > {std::make_pair("dedx_p_430.root", "430 mm"),
                                                          std::make_pair("dedx_p_530.root", "530 mm")});

  int c_arr[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  int m_arr[9] = {20, 20, 20,204, 20, 20, 20, 20, 20};
  int l_arr[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

  for (auto scan:list) {
    TString name = scan[0].first(0, 6);
    auto mg_r = new TMultiGraph();
    auto mg_c = new TMultiGraph();
    TLegend* l_r = new TLegend(gStyle->GetPadLeftMargin(),
                          1-gStyle->GetPadTopMargin() - 0.15,
                          gStyle->GetPadLeftMargin() + 0.3,
                          1-gStyle->GetPadTopMargin());
    TLegend* l_c = new TLegend(gStyle->GetPadLeftMargin(),
                          1-gStyle->GetPadTopMargin() - 0.15,
                          gStyle->GetPadLeftMargin() + 0.3,
                          1-gStyle->GetPadTopMargin());
    mg_r->SetName(name + "_resol");
    mg_c->SetName(name + "_charge");

    // scan iterator (for plotting)
    int i = 0;

    float max_y_r = 0;
    float min_y_r = 1e6;
    float max_y_c = 0;
    float min_y_c = 1e6;
    for (auto file:scan) {
      auto file_in = new TFile(input_prefix + file.first, "READ");
      auto graph_resol = (TGraphErrors*)file_in->Get("resol");
      auto graph_charg = (TGraphErrors*)file_in->Get("charge");
      graph_resol->SetMarkerColor(c_arr[i]);
      graph_charg->SetMarkerColor(c_arr[i]);
      graph_resol->SetLineColor(c_arr[i]);
      graph_charg->SetLineColor(c_arr[i]);
      graph_resol->SetMarkerStyle(m_arr[i]);
      graph_charg->SetMarkerStyle(m_arr[i]);
      graph_resol->SetLineStyle(l_arr[i]);
      graph_charg->SetLineStyle(l_arr[i]);
      mg_r->Add(graph_resol, "p");
      mg_c->Add(graph_charg, "p");

      l_r->AddEntry(graph_resol, file.second);
      l_c->AddEntry(graph_charg, file.second);


      if (TMath::MaxElement(graph_resol->GetN(), graph_resol->GetY()) > max_y_r) max_y_r = TMath::MaxElement(graph_resol->GetN(), graph_resol->GetY());
      if (TMath::MinElement(graph_resol->GetN(), graph_resol->GetY()) < min_y_r) min_y_r = TMath::MinElement(graph_resol->GetN(), graph_resol->GetY());
      if (TMath::MaxElement(graph_charg->GetN(), graph_charg->GetY()) > max_y_c) max_y_c = TMath::MaxElement(graph_charg->GetN(), graph_charg->GetY());
      if (TMath::MinElement(graph_charg->GetN(), graph_charg->GetY()) < min_y_c) min_y_c = TMath::MinElement(graph_charg->GetN(), graph_charg->GetY());
      mg_r->GetYaxis()->SetTitle(graph_resol->GetYaxis()->GetTitle());
      mg_c->GetYaxis()->SetTitle(graph_charg->GetYaxis()->GetTitle());
      mg_r->GetXaxis()->SetTitle(graph_resol->GetXaxis()->GetTitle());
      mg_c->GetXaxis()->SetTitle(graph_charg->GetXaxis()->GetTitle());
      ++i;
    }
    mg_r->SetMaximum(1.1*max_y_r);
    mg_r->SetMinimum(0.9*min_y_r);
    mg_c->SetMaximum(1.01*max_y_c);
    mg_c->SetMinimum(0.99*min_y_c);
    mg_c->GetXaxis()->ZoomOut(2);
    mg_r->GetXaxis()->ZoomOut(2);

    out_file->cd();
    mg_r->Write();
    mg_c->Write();

    c1.cd();
    mg_r->Draw("a");
    l_r->Draw();
    c1.Write("c_resol");
    mg_c->Draw("a");
    l_c->Draw();
    c1.Write("c_mean");
    //c1.WaitPrimitive();
  }

  out_file->Close();

  exit(0);
}