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

  //Initialize tree
  outtree = new TTree("outtree","outtree");
  outtree->Branch("ev",         &_ev,         "ev/I");
  outtree->Branch("dEdx",       &_dEdx,       "dEdx/F");
  outtree->Branch("npoints",    &_npoints,    "npoints/I");
  outtree->Branch("angle_xy",   &_angle_xy,   "_angle_xy/F");
  outtree->Branch("angle_yz",   &_angle_yz,   "_angle_yz/F");

  outtree->Branch("multiplicity",
                  &_multiplicity,
                  TString::Format("multiplicity[%i]/I", geom::nPadx)
                  );
  outtree->Branch("mult_rob",
                  &_multiplicity_robust,
                  TString::Format("multiplicity_robust[%i]/I", geom::nPadx)
                  );
  outtree->Branch("charge",
                  &_charge,
                  TString::Format("_charge[%i]/I", geom::nPadx)
                  );
  outtree->Branch("maxcharge_frac",
                  &_maxcharge_frac,
                  TString::Format("_maxcharge_frac[%i]/F", geom::nPadx)
                  );
  outtree->Branch("maxcharge_time",
                  &_maxcharge_time,
                  TString::Format("_maxcharge_time[%i]/I", geom::nPadx)
                  );

  outtree->Branch("pad_charge",
                  &_pad_charge,
                  TString::Format("_pad_charge[10][%i]/I", geom::nPadx)
                  );
  outtree->Branch("pad_time",
                  &_pad_time,
                  TString::Format("_pad_time[10][%i]/I", geom::nPadx)
                  );

  outtree->Branch("wf_width",
                  &_wf_width,
                  TString::Format("_wf_width[10][%i]/I", geom::nPadx)
                  );

  outtree->Branch("wf_fwhm",
                  &_wf_fwhm,
                  TString::Format("_wf_fwhm[10][%i]/I", geom::nPadx)
                  );

  outtree->Branch("pad_x",
                  &_pad_x,
                  TString::Format("_pad_x[10][%i]/I", geom::nPadx)
                  );
  outtree->Branch("pad_y",
                  &_pad_y,
                  TString::Format("_pad_y[10][%i]/I", geom::nPadx)
                  );

  _output_vector.push_back(outtree);

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

  _delta_t_fst = new TH1F("delta_t_fst", "", 300, -150., 150.);
  _delta_t_scd = new TH1F("delta_t_scd", "", 300, -150., 150.);

  _max_charge_time  = new TH1F("max_charge_time","",511,0,511);
  _max_charge_pos   = new TH1F("max_charge_pos", "", 40, 0., 40.);

  _XZ_leading = new TH2F("XZ_leading", "", geom::Nsamples, 0, geom::Nsamples, geom::nPadx, 0, geom::nPadx);
  _XZ_bias   = new TH1F("XZ_bias", "", 500, -0.05, 0.05);

  _angle = new TH2F("angle", "YZ angle vs XY angle", 150, 0., 1.5, 150, 0., 1.5);

  _delta_t_angle = new TH2F("dt_angle", "", 300, -150., 150., 150, 0., 1.5);

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

  _output_vector.push_back(_max_charge_time);
  _output_vector.push_back(_max_charge_pos);

  _output_vector.push_back(_XZ_leading);
  _output_vector.push_back(_XZ_bias);

  _selEvents = 0;

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  return true;
}

bool dEdxAna::ProcessEvent(const TEvent *event) {

  for(int trkID=0; trkID<(int)event->GetTracks().size(); trkID++){
    // add first digit with a track number
    _ev = event->GetID() + (trkID+1) * 1e8;
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
    // vector of charge in a cluster
    std::vector <double> QsegmentS;
    auto cols = itrack->GetCols(_invert);
    for(auto col:cols) if(col.size()){
      int colQ = 0;
      auto it_x = col[0]->GetCol(_invert);
      std::vector <THit*> Qpads;
      int z_max, x_max, q_max;
      z_max = x_max = q_max = 0;

      auto robust_col = GetRobustPadsInColumn(col);
      for(auto h:robust_col){
        colQ+=h->GetQ();
        _hTime->Fill(h->GetTime());
        Qpads.push_back(h);
        // study thw WF
        if (h->GetQ() > q_max) {
          q_max = h->GetQ();
          z_max = h->GetTime();
          x_max = h->GetCol(_invert);
        }
      } // loop over rows

      _maxcharge_frac[it_x] = (float) q_max/colQ;
      _maxcharge_time[it_x] = z_max;
      _charge[it_x] = colQ;

      if (colQ) {
        _XZ_leading->Fill(z_max, x_max);
        QsegmentS.push_back(colQ);
      }
      _un_trunk_cluster->Fill(colQ);

      _mult->Fill(col.size());

      _multiplicity[it_x] = col.size();
      _multiplicity_robust[it_x] = robust_col.size();

      sort(Qpads.begin(), Qpads.end(), [](THit* x1,
                                          THit* x2) {
                                            return x1->GetQ() > x2->GetQ();
                                          });

      for (uint pad_id = 0; pad_id < 10; pad_id++) {
        if (pad_id < Qpads.size()) {
          _pad_time[pad_id][it_x]   = Qpads[pad_id]->GetTime();
          _pad_charge[pad_id][it_x] = Qpads[pad_id]->GetQ();
          _wf_width[pad_id][it_x]   = Qpads[pad_id]->GetWidth();
          _wf_fwhm[pad_id][it_x]    = Qpads[pad_id]->GetFWHM();
          _pad_x[pad_id][it_x]      = Qpads[pad_id]->GetCol();
          _pad_y[pad_id][it_x]      = Qpads[pad_id]->GetRow();
        } else {
          _pad_time[pad_id][it_x]   = -9999;
          _pad_charge[pad_id][it_x] = -9999;
          _wf_width[pad_id][it_x]   = -9999;
          _wf_fwhm[pad_id][it_x]    = -9999;
          _pad_x[pad_id][it_x]      = -9999;
          _pad_y[pad_id][it_x]      = -9999;
        }
      }
    } // loop over column

    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    float CT= totQ / (alpha * QsegmentS.size());
    _hdEdx->Fill(CT);

    _npoints = QsegmentS.size();
    _dEdx = CT;

    outtree->Fill();

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
  } // loop over tracks
  return true;
}

bool dEdxAna::WriteOutput() {
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
