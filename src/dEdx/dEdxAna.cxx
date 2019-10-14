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

  _max_charge_pad = new TH1F("MaxChargePad", "", 1000, 0., 10000);

  _mult_graph = new TGraphErrors();
  _mult_graph->SetName("mult_graph");

  _mult_graph_err = new TGraphErrors();
  _mult_graph_err->SetName("mult_graph_err");

  _un_trunk_cluster = new TH1F("un_trunc_cluster", "", 1000, 0., 10000);

  _fst_pad_charge = new TH1F("fst_pad", "", 1000, 0., 10000);
  _scd_pad_charge = new TH1F("scd_pad", "", 1000, 0., 10000);
  _trd_pad_charge = new TH1F("trd_pad", "", 1000, 0., 10000);
  _fth_pad_charge = new TH1F("fth_pad", "", 1000, 0., 10000);

  _max_charge_time  = new TH1F("max_charge_time","",511,0,511);
  _max_charge_pos   = new TH1F("max_charge_pos", "", 40, 0., 40.);

  _XZ_leading = new TH2F("XZ_leading", "", geom::Nsamples, 0, geom::Nsamples, geom::nPadx, 0, geom::nPadx);

  for (auto i = 0; i < 4; ++i) {
    _charge_per_mult.push_back(new TH1F(Form("charge_per_mult_%i", i), "", 1000, 0., 10000.));
  }

  _output_vector.push_back(_hdEdx);
  _output_vector.push_back(_hTime);
  _output_vector.push_back(_mult);

  _output_vector.push_back(_mult_graph);
  _output_vector.push_back(_mult_graph_err);
  _output_vector.push_back(_max_charge_pad);
  _output_vector.push_back(_un_trunk_cluster);

  _output_vector.push_back(_fst_pad_charge);
  _output_vector.push_back(_scd_pad_charge);
  _output_vector.push_back(_trd_pad_charge);
  _output_vector.push_back(_fth_pad_charge);

  _output_vector.push_back(_max_charge_time);
  _output_vector.push_back(_max_charge_pos);

  _output_vector.push_back(_XZ_leading);

  for (uint i = 0; i < _charge_per_mult.size(); ++i)
    _output_vector.push_back(_charge_per_mult[i]);

  for (auto i = 0; i < geom::nPadx; ++i) {
    _mult_col[i] = new TH1F(Form("Mult_col_%i", i), "multiplicity", 10, 0., 10.);
    _output_vector.push_back(_mult_col[i]);
  }
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
    if(sel::GetColsMaxSep(itrack)>8) return false;
    if (sel::GetColsMaxGap(itrack) > 0) return false;
    std::vector<double> fit_v = sel::GetFitParams(itrack);
    if(fit_v[0]>1.0e6) return false;
    if (abs(fit_v[2]) > 0.08) return false;

    _store_event = true;

    //sel::Get3DFitParams(itrack);

    //If survives the selection, use track info:
    _selEvents++;

    if (_test_mode)
      if(_selEvents%10 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
    std::vector <double> QsegmentS;
    for(auto col:itrack->GetCols()) if(col.size()){
      int colQ = 0;
      auto it_x = col[0]->GetCol();
      std::vector <double> Qpads;
      int z_max, x_max, q_max;
      q_max = 0;
      for(auto h:col){
        colQ+=h->GetQ();
        _hTime->Fill(h->GetTime());
        Qpads.push_back(h->GetQ());
        if (h->GetQ() > q_max) {
          q_max = colQ;
          z_max = h->GetTime();
          x_max = h->GetCol();
        }
      }
      if (colQ) _XZ_leading->Fill(z_max, x_max);
      if (colQ) QsegmentS.push_back(colQ);
      _un_trunk_cluster->Fill(colQ);

      if (col.size() != 0 && col.size() <= _charge_per_mult.size())
        _charge_per_mult[col.size()-1]->Fill(colQ);

      _mult->Fill(col.size());
      _mult_col[it_x]->Fill(col.size());

      sort(Qpads.begin(), Qpads.end(), [](double x1, double x2){return x1 > x2;});
      if (Qpads.size() > 0) _fst_pad_charge->Fill(Qpads[0]);
      if (Qpads.size() > 1) _scd_pad_charge->Fill(Qpads[1]);
      if (Qpads.size() > 2) _trd_pad_charge->Fill(Qpads[2]);
      if (Qpads.size() > 3) _fth_pad_charge->Fill(Qpads[3]);
    } // loop over column
    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    double dEdx= totQ / (alpha * QsegmentS.size());
    _hdEdx->Fill(dEdx);

    // look for max charge in the pad in the event
    std::vector<THit*> hits = itrack->GetHits();
    auto it_max = std::max_element(hits.begin(), hits.end(),
                  [](THit* h1, THit* h2) { return h1->GetQ() < h2->GetQ(); });

    if (it_max != hits.end()) {
      int MaxCharge       = (*it_max)->GetQ();
      int MaxCharge_pos   = (*it_max)->GetCol();
      int MaxCharge_time  = (*it_max)->GetTime();
      _max_charge_pad->Fill(MaxCharge);
      _max_charge_pos->Fill(MaxCharge_pos);
      _max_charge_time->Fill(MaxCharge_time);
    }

    if(_batch == 0) {
      DrawCharge();
      DrawSelection(event,trkID);
    }
  }
  return true;
}

bool dEdxAna::WriteOutput() {
  std::cout << "Process result for output................";
  for (auto i = 0; i < geom::nPadx; ++i) {
    _mult_graph->SetPoint(_mult_graph->GetN(), i, _mult_col[i]->GetMean());
    _mult_graph->SetPointError(_mult_graph->GetN() - 1, 0, _mult_col[i]->GetRMS());

    _mult_graph_err->SetPoint(_mult_graph_err->GetN(), i, _mult_col[i]->GetMean());
    _mult_graph_err->SetPointError(_mult_graph_err->GetN() - 1, 0, _mult_col[i]->GetMeanError());
  }
  std::cout << "done" << std::endl;
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

bool dEdxAna::DrawCharge() {
  TCanvas* c1 = new TCanvas("c1", "c1", 0., 600., 800., 600.);
  _un_trunk_cluster->Draw();
  c1->Update();
  return true;
}

int main(int argc, char** argv) {
  auto ana = new dEdxAna(argc, argv);
  ana->Initialize();
  ana->Loop(ana->GetEventList());
  ana->WriteOutput();

  return 0;
}
