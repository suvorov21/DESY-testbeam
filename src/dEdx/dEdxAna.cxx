//#include <iostream>

#include "dEdxAna.hxx"

dEdxAna::dEdxAna(int argc, char** argv): AnalysisBase(argc, argv) {
  (void)argc;
}

bool dEdxAna::Initialize() {
  AnalysisBase::Initialize();

  _hdEdx  = new TH1F("dEdx","",300,0,5000);
  _hTime  = new TH1F("tDist","",511,0,511);
  _mult   = new TH1F("Mult", "multiplicity", 10, 0., 10.);
  _output_vector.push_back(_hdEdx);
  _output_vector.push_back(_hTime);
  _output_vector.push_back(_mult);
  _selEvents = 0;

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  return true;
}

bool dEdxAna::ProcessEvent(const TEvent *event) {
  double alpha = 0.625;
  for(int trkID=0; trkID<(int)event->GetTracks().size(); trkID++){
    TTrack* itrack = event->GetTracks()[trkID];
    if(_verbose == 2){
      std::cout << "sel::GetNonZeroCols(event,trkID).size(): " << sel::GetNonZeroCols(itrack).size() << std::endl;
      std::cout << "sel::GetNonZeroRows(event,trkID).size(): " << sel::GetNonZeroRows(itrack).size() << std::endl;
    }
    if(sel::GetNonZeroCols(itrack).size() != 36) return false;
    if(sel::GetColsMaxSep(itrack)>5) return false;
    if(sel::GetFitParams(itrack)[0]>1.0e6) return false;

    //sel::Get3DFitParams(itrack);

    //If survives the selection, use track info:
    _selEvents++;
    if(_batch == 0) DrawSelection(event,trkID);
    if (_test_mode)
      if(_selEvents%10 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
    std::vector <double> QsegmentS;
    for(auto col:itrack->GetCols()) if(col.size()){
      int colQ = 0;
      for(auto h:col){
        colQ+=h->GetQ();
        _hTime->Fill(h->GetTime());
      }
      if(colQ) QsegmentS.push_back(colQ);
      _mult->Fill(col.size());
    }
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
  AnalysisBase::WriteOutput();

  // std::cout << "selEvents: " << _selEvents << std::endl;
  // std::cout << "Write dedx output........................";

  // if (!_file_out)
  //   return true;

  // auto file = new TFile(_file_out_name.Data(), "UPDATE");
  // // write
  // file->Close();

  // std::cout << "done" << std::endl;
  return true;
}

int main(int argc, char** argv) {
  auto ana = new dEdxAna(argc, argv);
  ana->Initialize();
  ana->Loop(ana->GetEventList());
  ana->WriteOutput();

  return 0;
}
