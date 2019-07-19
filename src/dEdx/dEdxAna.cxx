//#include <iostream>

#include "dEdxAna.hxx"

dEdxAna::dEdxAna(int argc, char** argv): AnalysisBase(argc, argv) {
  (void)argc;
}

bool dEdxAna::Initialize() {
  AnalysisBase::Initialize();

  _hdEdx = new TH1F("dEdx","",300,0,5000);
  _output_vector.push_back(_hdEdx);
  _selEvents = 0;
  _verbose = 0;

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  return true;
}

bool dEdxAna::ProcessEvent(const Event event) {
  double alpha = 0.625;
  for(int trkID=0; trkID<event.tracks.size(); trkID++){
    if(sel::GetNonZeroCols(event,trkID).size() != 36) return false;
    if(sel::GetNonZeroRows(event,trkID).size()>5) return false;
    if(sel::GetFitQuality(event,trkID)>1.0e6) return false;
    //If survives the selection, use track info:
    _selEvents++;
    DrawSelection(event,trkID);
    if(_selEvents%10 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
    std::vector <double> QsegmentS =  sel::GetNonZeroCols(event,trkID);   
    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    double dEdx= totQ / (alpha * QsegmentS.size());
    _hdEdx->Fill(dEdx);
  }
  return true;
}

bool dEdxAna::WriteOutput() {
  std::cout << "Write spatial output.....................";
  AnalysisBase::WriteOutput();

  std::cout << "selEvents: " << _selEvents << std::endl;

  if (!_file_out)
    return true;

  auto file = new TFile(_file_out_name.Data(), "UPDATE");
  // write
  file->Close();

  std::cout << "done" << std::endl;
  return true;
}

int main(int argc, char** argv) {
  auto ana = new dEdxAna(argc, argv);
  ana->Initialize();
  ana->Loop(ana->GetEventList());
  ana->WriteOutput();

  return 1;
}
