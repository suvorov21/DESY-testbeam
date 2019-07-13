// ROOT
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

// c++
#include <iostream> // stream
#include <unistd.h> // getopt on Mac
#include <fstream>  // read file lists

// project
#include "../utils/SetT2KStyle.hxx"
#include "../utils/Geom.hxx"

void help(std::string name);

int main(int argc, char** argv) {
  // ini params
  auto batch        = false;
  auto test_mode    = false;
  auto verbose      = 1;

  TString file_in_name  = "";
  TString file_out_name = "";

  // read CLI
  for (;;) {
    int c = getopt(argc, argv, "i:o:bvdm");
    if (c < 0) break;
    switch (c) {
      case 'i' : file_in_name     = optarg;       break;
      case 'b' : batch            = true;         break;
      case 'v' : verbose          = atoi(optarg); break;
      case 'd' : test_mode        = true;         break;
      case 'm' : help(argv[0]);                   break;
      case '?' : help(argv[0]);
    }
  }
  if (verbose) {
    std::cout << "Running      " << argv[0] << std::endl;
    std::cout << "Verbosity  : " << verbose << std::endl;
    std::cout << "Test mode  : " << test_mode << std::endl;
    std::cout << "Batch mode : " << batch << std::endl;
  }

  // define style
  auto T2KstyleIndex = 2;
  // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
  auto localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper

  auto t2kstyle = SetT2KStyle(T2KstyleIndex, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());
  gROOT->ForceStyle();

  TApplication* app;
  if (!batch)
    app = new TApplication("App", &argc, argv);

  auto c1         = new TCanvas("c1", "2D max ampl",  0, 0, plot::canvas_x_size, plot::canvas_y_size);
  auto c2         = new TCanvas("c2", "3D",           plot::canvas_x_size, 0, plot::canvas_x_size, plot::canvas_y_size);
  auto h2d        = new TH2F("h2d",   "",
    geom::nPadx+2, -1, geom::nPadx + 1, geom::nPady+2, -1, geom::nPady+1);
  auto h3d        = new TH3F("h3d",   "",
    geom::nPadx+2, -1, geom::nPadx + 1, geom::nPady+2, -1, geom::nPady+1, geom::Nsamples+2, -1, geom::Nsamples+1);

  // read data
  Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];

  auto file_in = new TFile(file_in_name.Data(), "READ");
  if (!file_in->IsOpen()) {
    std::cerr << "Input file is not open. " << file_in_name << std::endl;;
    exit(1);
  }
  // read and chain files
  auto chain = new TChain("chain");

  if (file_in_name.Contains(".root")) {
    std::cout << "adding filename" <<" " << file_in_name << std::endl;
    chain->AddFile(file_in_name);
  } else {
    std::ifstream fList(file_in_name.Data());
    if (!fList.good()) {
      std::cerr << "Can not read input " << file_in_name << std::endl;
      exit(1);
    }
    while (fList.good()) {
      std::string filename;
      getline(fList, filename);
      if (fList.eof()) break;
      chain->AddFile(filename.c_str());
    }
  }


  chain->SetBranchAddress("PadAmpl", padAmpl);
  Int_t N_events = chain->GetEntries();
  if (test_mode)
    N_events = std::min(N_events, 30);

  if (verbose == 1)
    std::cout << "Processing" << std::endl;

  for (auto eventID = 0; eventID < N_events; ++eventID) {
    chain->GetEntry(eventID);

    std::cout << "Event " << eventID << std::endl;

    h3d->Reset();
    h2d->Reset();

    // SELECTIONS IF NECESSARY

    for (auto x = 0; x < geom::nPadx; ++x) {
      for (auto y = 0; y < geom::nPady; ++y) {
        auto maxAmple = 0;

        for (auto it = 0; it < geom::Nsamples; ++it) {
          if (padAmpl[x][y][it])
            h3d->Fill(x, y, it, padAmpl[x][y][it]);
          if (padAmpl[x][y][it] > maxAmple)
            maxAmple = padAmpl[x][y][it];
        } // end of loop over time

        if (maxAmple == 0)
          continue;

        h2d->Fill(x, y, maxAmple);

        // PUT YOUR ANALYSIS HERE

      } // end of loop over y
    } // end of loop over x
    if (!batch) {
      c1->cd();
      h2d->Draw("colz");
      c2->cd();
      h3d->Draw("box1");
      c2->WaitPrimitive();
    }
  } // end of loop over events

  // Draw in case of non-batch run
  if (!batch) {
    exit(1);
  }

  std::cout << "Complete" << std::endl;
  if (!batch)
    app->Run();

  return 1;
}

void help(std::string name)
{
  std::cout << name << " usage\n" << std::endl;
  std::cout << "   -i <input_file>      : input file name with a path" << std::endl;
  std::cout << "   -b                   : run in batch mode" << std::endl;
  std::cout << "   -v <verbose_lvel>    : verbosity level" << std::endl;
  std::cout << "   -d                   : test mode. run over first 30 events" << std::endl;
  std::cout << "   -h                   : print ROOT help" << std::endl;
  std::cout << "   -m                   : print " << name << " help" << std::endl;
  exit(1);
}