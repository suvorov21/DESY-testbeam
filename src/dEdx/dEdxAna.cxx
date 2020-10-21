//#include <iostream>

#include "dEdxAna.hxx"

dEdxAna::dEdxAna(int argc, char** argv): AnalysisBase(argc, argv) {
  (void)argc;
}

bool dEdxAna::Initialize() {
  AnalysisBase::Initialize();

  /**
   * Tree declaration
   */

  /// tree with entry per cluster
  _cluster_tree = new TTree("cluster_tree", "cluster_tree");
  _cluster_tree->Branch("event_id",    &_event_id);
  _cluster_tree->Branch("mult",        &_multiplicity);
  _cluster_tree->Branch("dEdx_trunc",  &_dedx_truncated);
  _cluster_tree->Branch("angle_xy",    &_angle_xy);
  _cluster_tree->Branch("angle_yz",    &_angle_yz);
  _cluster_tree->Branch("charge",      &_charge,    "_charge[10]/I");
  _cluster_tree->Branch("time",        &_time,      "_time[10]/I");

  _output_vector.push_back(_cluster_tree);


  /**
   * Histogram declaration
   */

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

  _delta_t_fst = new TH1F("delta_t_fst", "", 300, -150., 150.);
  _delta_t_scd = new TH1F("delta_t_scd", "", 300, -150., 150.);

  _max_charge_time  = new TH1F("max_charge_time","",511,0,511);
  _max_charge_pos   = new TH1F("max_charge_pos", "", 40, 0., 40.);

  _XZ_leading = new TH2F("XZ_leading", "", geom::Nsamples, 0, geom::Nsamples, geom::nPadx, 0, geom::nPadx);
  _XZ_bias   = new TH1F("XZ_bias", "", 500, -0.05, 0.05);

  _angle = new TH2F("angle", "YZ angle vs XY angle", 150, 0., 1.5, 150, 0., 1.5);

  _delta_t_angle = new TH2F("dt_angle", "", 300, -150., 150., 150, 0., 1.5);

  for (auto i = 0; i < 4; ++i) {
    _charge_per_mult.push_back(new TH1F(Form("charge_per_mult_%i", i), "", 1000, 0., 10000.));
  }

  _output_vector.push_back(_hdEdx);
  _output_vector.push_back(_hTime);
  _output_vector.push_back(_mult);

  _output_vector.push_back(_angle);

  _output_vector.push_back(_delta_t_fst);
  _output_vector.push_back(_delta_t_scd);
  _output_vector.push_back(_delta_t_angle);

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
  _output_vector.push_back(_XZ_bias);

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

  for(int trkID=0; trkID<(int)event->GetTracks().size(); trkID++){
    // add first digit with a track number
    _event_id = event->GetID() + (trkID+1) * 1e8;
    TTrack* itrack = event->GetTracks()[trkID];
    if(_verbose > 1){
      std::cout << "sel::GetNonZeroCols(event,trkID).size(): ";
      std::cout << sel::GetNonZeroCols(itrack, _invert).size() << std::endl;
      std::cout << "sel::GetColsMaxSep(event,trkID).size():  ";
      std::cout << sel::GetColsMaxSep(itrack, _invert) << std::endl;
      std::cout << "sel::GetColsMaxGap(event,trkID).size():  ";
      std::cout << sel::GetColsMaxGap(itrack, _invert) << std::endl;
    }

    if (!sel::CrossingTrackSelection(itrack, _invert, _verbose))
      continue;
    std::vector<double> fit_v = sel::GetFitParams(itrack, _invert);
    std::vector<double> fit_xz = sel::GetFitParamsXZ(itrack, _invert);

    _angle->Fill(abs(fit_v[2]), abs(fit_xz[2] * sel::v_drift_est));
    _angle_xy = fit_v[2];
    _angle_yz = fit_xz[2] * sel::v_drift_est;

    _store_event = true;

    // if survives the selection, use track info:
    _selEvents++;

    if (_test_mode)
      if(_selEvents%10 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
    std::vector <double> QsegmentS;
    for(auto col:itrack->GetCols(_invert)) if(col.size()){
      int colQ = 0;
      auto it_x = col[0]->GetCol(_invert);
      std::vector <std::pair<int, int> > Qpads;
      int z_max, x_max, q_max;
      z_max = x_max = q_max = 0;
      for(auto h:col){
        colQ+=h->GetQ();
        _hTime->Fill(h->GetTime());
        Qpads.push_back(std::make_pair(h->GetTime(), h->GetQ()));
        if (h->GetQ() > q_max) {
          q_max = h->GetQ();
          z_max = h->GetTime();
          x_max = h->GetCol(_invert);
        }
      }
      if (colQ) {
        _XZ_leading->Fill(z_max, x_max);
        QsegmentS.push_back(colQ);
      }
      _un_trunk_cluster->Fill(colQ);

      if (col.size() != 0 && col.size() <= _charge_per_mult.size())
        _charge_per_mult[col.size()-1]->Fill(colQ);

      _mult->Fill(col.size());
      _mult_col[it_x]->Fill(col.size());

      sort(Qpads.begin(), Qpads.end(), [](std::pair<int, int> x1,
                                          std::pair<int, int> x2) {
                                            return x1.second > x2.second;
                                          });

      if (Qpads.size() > 0) _fst_pad_charge->Fill(Qpads[0].second);
      if (Qpads.size() > 1) {
        _scd_pad_charge->Fill(Qpads[1].second);
        _delta_t_fst->Fill(Qpads[1].first - Qpads[0].first);
        _delta_t_angle->Fill(Qpads[1].first - Qpads[0].first,
                             abs(fit_xz[2] * sel::v_drift_est)
                             );
      }
      if (Qpads.size() > 2) {
        _trd_pad_charge->Fill(Qpads[2].second);
        _delta_t_scd->Fill(Qpads[2].first - Qpads[0].first);
      }
      if (Qpads.size() > 3) _fth_pad_charge->Fill(Qpads[3].second);

      Qpads.resize(10, std::make_pair(0, 0));
      for (auto padID = 0; padID < 10; ++padID) {
        _charge[padID] = Qpads[padID].second;
        _time[padID] = Qpads[padID].first;
      }

      _multiplicity = col.size();
      _cluster_tree->Fill();
    } // loop over column

    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    double dEdx= totQ / (alpha * QsegmentS.size());
    _hdEdx->Fill(dEdx);
    _dedx_truncated = dEdx;

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

    // if(_batch == 0) {
    //   DrawCharge();
    //   DrawSelection(event,trkID);
    // }
  } // loop over tracks
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
