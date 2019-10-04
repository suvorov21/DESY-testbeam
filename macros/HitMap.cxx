// ROOT
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
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
      case 'o' : file_out_name    = optarg;       break;
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
    std::cout << std::endl;
    std::cout << "Input  : " << file_in_name << std::endl;
    std::cout << "Output : " << file_out_name << std::endl;
  }

  // define style
  auto T2KstyleIndex = 2;
  // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
  auto localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper

  auto t2kstyle = T2K().SetT2KStyle(T2KstyleIndex, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());
  gROOT->ForceStyle();

  TApplication* app;
  if (!batch)
    app = new TApplication("App", &argc, argv);

  // define output
  file_out_name += "/HitMap.root";
  auto file_out   = new TFile(file_out_name.Data(), "RECREATE");
  if (!file_out->IsOpen()) {
    std::cerr << "Output file is not open. " << file_out_name << std::endl;;
    exit(1);
  }
  auto c1         = new TCanvas("c1", "", plot::canvas_x_size, plot::canvas_y_size);
  auto Nhits_h    = new TH2F("Nhits",   "", geom::nPadx+2, -1, geom::nPadx + 1, geom::nPady+2, -1, geom::nPady+1);
  auto MaxAmpl_h  = new TH2F("MaxAmpl", "", geom::nPadx+2, -1, geom::nPadx + 1, geom::nPady+2, -1, geom::nPady+1);

  auto charge     = new TH1F("charge", "", 4500, 0., 4500);

  // read data
  Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];

  // read and chain files
  auto chain = new TChain("tree");

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

  if (verbose == 1) {
    std::cout << "Processing" << std::endl;
    std::cout << "[                    ]   Nevents = " << N_events << "\r[";
  }
  for (auto eventID = 0; eventID < N_events; ++eventID) {
    chain->GetEntry(eventID);

    if (verbose == 1 && (eventID%(N_events/20)) == 0)
      std::cout << "." << std::flush;

    for (auto x = 0; x < geom::nPadx; ++x) {
      for (auto y = 0; y < geom::nPady; ++y) {
        auto maxAmple = 0;

        for (auto it = 0; it < geom::Nsamples; ++it) {
          ///if (it < 100 || it > 150)
          //  continue;
          if (padAmpl[x][y][it] > maxAmple)
            maxAmple = padAmpl[x][y][it];
        } // end of loop over time

        if (maxAmple == 0)
          continue;

        charge->Fill(maxAmple);

        Nhits_h->Fill(x, y);
        MaxAmpl_h->Fill(x, y, maxAmple);

      } // end of loop over y
    } // end of loop over x
  } // end of loop over events

  if (verbose == 1)
    std::cout << "]" << std::endl;

  auto AvAmpl_h = (TH2F*)MaxAmpl_h->Clone("AvgAmpl");
  AvAmpl_h->Divide(Nhits_h);

  file_out->cd();
  MaxAmpl_h->Write();
  Nhits_h->Write();
  AvAmpl_h->Write();

  charge->Write();

  file_out->Close();

  // Draw in case of non-batch run
  if (!batch) {
    c1->WaitPrimitive();
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
  std::cout << "   -o <output_path>     : output files path" << std::endl;
  std::cout << "   -b                   : run in batch mode" << std::endl;
  std::cout << "   -v <verbose_lvel>    : verbosity level" << std::endl;
  std::cout << "   -d                   : test mode. run over first 30 events" << std::endl;
  std::cout << "   -h                   : print ROOT help" << std::endl;
  std::cout << "   -m                   : print " << name << " help" << std::endl;
  exit(1);
}