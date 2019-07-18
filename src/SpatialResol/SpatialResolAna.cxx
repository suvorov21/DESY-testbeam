//#include <iostream>

#include "SpatialResolAna.hxx"

SpatialResolAna::SpatialResolAna(int argc, char** argv): AnalysisBase(argc, argv) {
  _iteration = -1;

  for (;;) {
    int c = getopt(argc, argv, "i:o:bvdmt");
    if (c < 0) break;
    switch (c) {
      case 't' : _iteration     = atoi(optarg);    break;
    }
  }

  if (_iteration == -1) {
    std::cerr << "ERROR. SpatialResolAna::SpatialResolAna(). Iteration should be defined as a input param" << std::endl;
    exit(1);
  }

}

bool SpatialResolAna::Initialize() {
  AnalysisBase::Initialize();

  // Initilise selection
  _selection = new DBSCANSelection();
  _selection->Initialize();

  return true;
}

bool SpatialResolAna::ProcessEvent(const Event event) {

  return true;
}

bool SpatialResolAna::WriteOutput() {
  std::cout << "Write spatial output.....................";
  AnalysisBase::WriteOutput();

  if (!_file_out)
    return true;

  auto file = new TFile(_file_out_name.Data(), "UPDATE");
  // write
  file->Close();

  std::cout << "done" << std::endl;
  return true;
}

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna(argc, argv);
  ana->Initialize();
  ana->Loop(ana->GetEventList());
  ana->WriteOutput();

  return 1;
}
