#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"

// TODO move it finally to the utils
float GetAverage(TGraphErrors* h, float& RMS, float& mean_e);

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "ERROR! SR_test::main(): no input " << std::endl;
    exit(1);
  }

  int limit = atoi(argv[2]);

  std::cout << "SR_test. Open file " << argv[1] << std::endl;
  auto f = new TFile(argv[1], "READ");
  auto t = (TTree*)f->Get("outtree");
  auto resol_total = new TH1F("resol", "", 200, -0.004, 0.004);
  t->Project("resol_total", "residual");

  if (!resol_total || !resol_total->GetEntries()) {
    std::cerr << "ERROR! SR_test::main(): no entries in resol_total" << std::endl;
    exit(1);
  }

  resol_total->Fit("gaus");
  if (!resol_total->GetFunction("gaus")) {
    std::cerr << "ERROR! SR_test::main(): resol_total fit fail" << std::endl;
    exit(1);
  }

  auto sigma = resol_total->GetFunction("gaus")->GetParameter(2);
  if (1.e6 * sigma > limit) {
    std::cerr << "ERROR! SR_test::main(): too large spatial resolution"
    << "\t" << 1.e6 * sigma << std::endl;
    exit(1);
  }

  // TGraphErrors* resol_final = (TGraphErrors*)f->Get("resol_final");

  // if (!resol_final->GetN()) {
  //   std::cerr << "ERROR! SR_test::main(): no entries in resol_final" << std::endl;
  //   exit(1);
  // }

  // float RMS, mean_e;
  // float mean = GetAverage(resol_final, RMS, mean_e);

  // if (1.e6 * mean > limit) {
  //   std::cerr << "ERROR! SR_test::main(): too large spatial resolution"
  //   << "\t" << 1.e6 * mean << std::endl;
  //   exit(1);
  // }

  // if (1.e6 * RMS > 150) {
  //   std::cerr << "ERROR! SR_test::main(): too large SR RMS"
  //   << "\t" << 1.e6 * RMS << std::endl;
  //   exit(1);
  // }

  std::cout << "SR_test::main(): test passed" << std::endl;
  return 0;
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