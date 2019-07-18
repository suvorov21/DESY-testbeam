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

  // Initilise selection
  _selection = new DBSCANReconstruction();
  _selection->Initialize();

  return true;
}

bool dEdxAna::ProcessEvent(const Event event) {
  double alpha = 0.625;
  _selEvents++;
  if(_selEvents%100 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
  for(int num=0; num<event.trackNum; num++){
    std::vector <double> QsegmentS;   
    for (int itx = 0; itx < geom::nPadx; ++itx){
      double Qsegment = 0;      
      for (uint ity = 0; ity < geom::nPady; ++ity){
        if (!event.twoD[num][itx][ity]) continue;  
        Qsegment+=event.twoD[num][itx][ity];
      }
      if(Qsegment) QsegmentS.push_back(Qsegment);
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
