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

  //Initialize tree
  outtree = new TTree("outtree","outtree");
  outtree->Branch("ev",&ev,"ev/I");
  outtree->Branch("dEdx",&dEdx,"dEdx/F");
  outtree->Branch("npoints",&npoints,"npoints/I");
  outtree->Branch("angle_xy",&angle_xy,"angle_xy/F");
  outtree->Branch("angle_yz",&angle_yz,"angle_yz/F");
  outtree->Branch("delta_t_fst",&delta_t_fst,TString::Format("delta_t_fst[%i]/F", geom::nPadx));
  outtree->Branch("delta_t_scd",&delta_t_scd,TString::Format("delta_t_scd[%i]/F", geom::nPadx));
  outtree->Branch("multiplicity",&multiplicity,TString::Format("multiplicity[%i]/I", geom::nPadx));

  outtree->Branch("charge",&charge,TString::Format("charge[%i]/F", geom::nPadx));
  outtree->Branch("maxcharge_frac",&maxcharge_frac,TString::Format("maxcharge_frac[%i]/F", geom::nPadx));
  outtree->Branch("maxcharge_time",&maxcharge_time,TString::Format("maxcharge_time[%i]/F", geom::nPadx));

  outtree->Branch("fst_pad_charge",&fst_pad_charge,TString::Format("fst_pad_charge[%i]/F", geom::nPadx));
  outtree->Branch("scd_pad_charge",&scd_pad_charge,TString::Format("scd_pad_charge[%i]/F", geom::nPadx));
  outtree->Branch("trd_pad_charge",&trd_pad_charge,TString::Format("trd_pad_charge[%i]/F", geom::nPadx));
  outtree->Branch("fth_pad_charge",&fth_pad_charge,TString::Format("fth_pad_charge[%i]/F", geom::nPadx));

  outtree->Branch("fst_pad_time",&fst_pad_time,TString::Format("fst_pad_time[%i]/F", geom::nPadx));
  outtree->Branch("scd_pad_time",&scd_pad_time,TString::Format("scd_pad_time[%i]/F", geom::nPadx));
  outtree->Branch("trd_pad_time",&trd_pad_time,TString::Format("trd_pad_time[%i]/F", geom::nPadx));
  outtree->Branch("fth_pad_time",&fth_pad_time,TString::Format("fth_pad_time[%i]/F", geom::nPadx));




  _output_vector.push_back(outtree);

  // Initilise selection
  _reconstruction = new DBSCANReconstruction();
  _reconstruction->Initialize();

  return true;
}

bool dEdxAna::ProcessEvent(const TEvent *event) {
  double alpha = 0.625;
  for(int trkID=0; trkID<(int)event->GetTracks().size(); trkID++){
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

    // if(fit_v[0]>1.0e6) return false;

    _angle->Fill(abs(fit_v[2]), abs(fit_xz[2] * sel::v_drift_est));

    angle_xy = abs(fit_v[2]);
    angle_yz = abs(fit_xz[2]*sel::v_drift_est);

    _store_event = true;

    // cut f0r slope in XZ
    // if (_invert) {
    //   _XZ_bias->Fill(fit_xz[2]*0.0028*32);
    // } else
    //   _XZ_bias->Fill(fit_xz[2]*0.0028*geom::GetMaxColumn(_invert));

    //sel::Get3DFitParams(itrack);

    //If survives the selection, use track info:
    _selEvents++;

    if (_test_mode)
      if(_selEvents%10 == 0) std::cout << "selEvents: " << _selEvents << std::endl;
    std::vector <double> QsegmentS;
    for(auto col:itrack->GetCols(_invert)) if(col.size()){
      int colQ = 0;
      auto it_x = col[0]->GetCol(_invert);
      std::vector <std::pair<int, double> > Qpads;
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

      std::cout<<col.size()<<" "<<Qpads.size()<<std::endl;

      maxcharge_frac[it_x] = (float) q_max/colQ;
      maxcharge_time[it_x] = z_max;
      charge[it_x] = colQ;

      if (colQ) {
        _XZ_leading->Fill(z_max, x_max);
        QsegmentS.push_back(colQ);
      }
      _un_trunk_cluster->Fill(colQ);

      if (col.size() != 0 && col.size() <= _charge_per_mult.size())
        _charge_per_mult[col.size()-1]->Fill(colQ);

      _mult->Fill(col.size());
      _mult_col[it_x]->Fill(col.size());

      multiplicity[it_x] = col.size();

      sort(Qpads.begin(), Qpads.end(), [](std::pair<int, double> x1,
                                          std::pair<int, double> x2) {
                                            return x1.second > x2.second;
                                          });

      if (Qpads.size() == 1) {
	_fst_pad_charge->Fill(Qpads[0].second);

	fst_pad_charge[it_x] = Qpads[0].second;
	fst_pad_time[it_x] = Qpads[0].first;
	scd_pad_charge[it_x] = -9999;
	scd_pad_time[it_x] = -9999;
	trd_pad_charge[it_x] = -9999;
	trd_pad_time[it_x] = -9999;
      }
      else if (Qpads.size() == 2) {
        _scd_pad_charge->Fill(Qpads[1].second);
        _delta_t_fst->Fill(Qpads[1].first - Qpads[0].first);
        _delta_t_angle->Fill(Qpads[1].first - Qpads[0].first,
                             abs(fit_xz[2] * sel::v_drift_est)
                             );
	delta_t_fst[it_x] = Qpads[1].first - Qpads[0].first;
	delta_t_scd[it_x] = -99999;
	fst_pad_charge[it_x] = Qpads[0].second;
	fst_pad_time[it_x] = Qpads[0].first;
	scd_pad_charge[it_x] = Qpads[1].second;
	scd_pad_time[it_x] = Qpads[1].first;
	trd_pad_charge[it_x] = -9999;
	trd_pad_time[it_x] = -9999;
      }
      if (Qpads.size() == 3) {
        _trd_pad_charge->Fill(Qpads[2].second);
        _delta_t_scd->Fill(Qpads[2].first - Qpads[0].first);
	delta_t_fst[it_x] = Qpads[1].first - Qpads[0].first;
	delta_t_scd[it_x] = Qpads[2].first - Qpads[0].first;
	fst_pad_charge[it_x] = Qpads[0].second;
	fst_pad_time[it_x] = Qpads[0].first;
	scd_pad_charge[it_x] = Qpads[1].second;
	scd_pad_time[it_x] = Qpads[1].first;
	trd_pad_charge[it_x] = Qpads[2].second;
	trd_pad_time[it_x] = Qpads[2].first;
      }
      if (Qpads.size() > 3) _fth_pad_charge->Fill(Qpads[3].second);
    } // loop over column
    sort(QsegmentS.begin(), QsegmentS.end());
    double totQ = 0.;
    Int_t i_max = round(alpha * QsegmentS.size());
    for (int i = 0; i < std::min(i_max, int(QsegmentS.size())); ++i) totQ += QsegmentS[i];
    float CT= totQ / (alpha * QsegmentS.size());
    _hdEdx->Fill(CT);
    
    npoints = QsegmentS.size();
    dEdx = CT;
    ev = 1;
    
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

    // if(_batch == 0) {
    //   DrawCharge();
    //   DrawSelection(event,trkID);
    // }
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
