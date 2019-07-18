//#include <iostream>

#include "SpatialResolAna.hxx"

SpatialResolAna::SpatialResolAna(int argc, char** argv): AnalysisBase(argc, argv) {
  (void)argc;
}

bool SpatialResolAna::Initialize() {
  AnalysisBase::Initialize();

  // Initilise selection
  _reconstruction = new CrossingReconstruction();
  _reconstruction->Initialize();

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
