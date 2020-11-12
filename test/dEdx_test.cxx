#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "ERROR! dEdx_test::main(): no input " << std::endl;
    exit(1);
  }

  std::cout << "dEdx_test. Open file " << argv[1] << std::endl;
  auto f = new TFile(argv[1], "READ");
  auto h = (TH1F*)f->Get("dEdx");

  if (!h->GetEntries()) {
    std::cerr << "ERROR! dEdx_test::main(): no entries in dEdx " << std::endl;
    exit(1);
  }

  h->Fit("gaus");

  auto max    = h->GetMaximum();
  auto max_x  = h->GetBinCenter(h->GetMaximumBin());
  auto FWHM   = h->GetBinCenter(h->FindLastBinAbove(max/2))
              - h->GetBinCenter(h->FindFirstBinAbove(max/2));

  h->Fit("gaus", "", "", max_x-3*FWHM, max_x + 3*FWHM);
  TF1* fit = (TF1*)h->GetFunction("gaus");
  if (!fit) {
    std::cerr << "ERROR! dEdx_test::main(): no dEdx fit" << std::endl;
    exit(1);
  }

  auto mean     = fit->GetParameter(1);
  auto sigma    = fit->GetParameter(2);

  if (sigma / mean > 0.2) {
    std::cerr << "ERROR! dEdx_test::main(): large dEdx resolution" << std::endl;
    exit(1);
  }

  std::cout << "dEdx_test::main(): test passed" << std::endl;
  return 0;
}