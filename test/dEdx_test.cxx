#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "ERROR! dEdx_test::main(): no input " << std::endl;
    exit(1);
  }

  std::cout << "dEdx_test. Open file " << argv[1] << std::endl;
  auto f = TFile::Open(argv[1], "READ");
  auto t = (TTree*)f->Get("outtree");
  auto h = new TH1F("h", "", 200, 0., 4000.);
  t->Project("h", "dEdx", "dEdx > 0");

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
