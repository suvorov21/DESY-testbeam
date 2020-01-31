#include "TGraphErrors.h"
#include "TH2F.h"
#include "TH1F.h"


TGraphErrors* prf(TH2F* _PRF_histo) {
  TGraphErrors* _PRF_graph = new TGraphErrors();
  TAxis* ax = _PRF_histo->GetXaxis();

  TH1F* h = new TH1F("h", "", ax->GetNbins(), ax->GetBinLowEdge(1), ax->GetBinLowEdge(ax->GetNbins()) + ax->GetBinWidth(ax->GetNbins()));

  for (auto i = 1; i < _PRF_histo->GetXaxis()->GetNbins(); ++i) {

    TH1D* temp_h = _PRF_histo-> ProjectionY(Form("projections_bin_%i", i), i, i);

    double x = _PRF_histo->GetXaxis()->GetBinCenter(i);
    double y = temp_h->GetBinCenter(temp_h->GetMaximumBin());

    float start = -1.;
    float end   = -1.;
    float max = temp_h->GetMaximum();

    for (Int_t bin = 0; bin < temp_h->GetXaxis()->GetNbins(); ++bin) {
      if (start == -1. && temp_h->GetBinContent(bin) >= max / 2.)
        start = temp_h->GetBinCenter(bin);

      if (end == -1. && start != -1. && temp_h->GetBinContent(bin) <= max / 2.)
        end = temp_h->GetBinCenter(bin);
    }

    float e = end - start;

    _PRF_graph->SetPoint(_PRF_graph->GetN(), x, y);
    _PRF_graph->SetPointError(_PRF_graph->GetN()-1, 0, e/2.);

    if (abs(x) > 0.01)
      continue;

    h->SetBinContent(i, y);
  } // end of PRF histo profiling

  Float_t mean  = h->GetMean();
  Float_t sigma = 0.;
  float max   = h->GetMaximum();
  float start = h->GetBinLowEdge(h->FindFirstBinAbove(max/2));
  float end   = h->GetBinLowEdge(h->FindLastBinAbove(max/2)) +
  h->GetBinWidth(h->FindLastBinAbove(max/2));

  sigma = 0.5 * (end - start);

  std::cout << sigma*1e2 << std::endl;

  //TCanvas c1("c1", "2D max ampl",  0, 0, 800, 600);
  h->Draw();
  //c1.WaitPrimitive();


  return NULL;
}