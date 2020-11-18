#include "Math/Functor.h"
#include "Fit/Fitter.h"

#include "InclinedTracks.hxx"

//********************************************************************
InclinedTracks::InclinedTracks(int argc, char** argv):
//********************************************************************
  SpatialResolAna(argc, argv) {
    std::cout << "hey there" << std::endl;

      _circle_function_up = new TF1("circle_up",
    "-sqrt([0]*[0] - TMath::Power(x+0.198 - [1] * [0], 2)) + [0] * sqrt(1-[1]*[1]) + [2]",
    -0.5, 0.5);
  _circle_function_up->SetParName(0, "radius");
  _circle_function_up->SetParName(1, "sin(alpha)");
  _circle_function_up->SetParName(2, "target");

  _circle_function_dn = new TF1("circle_dn",
    "sqrt([0]*[0] - TMath::Power(x+0.198 - [1] * [0], 2)) - [0] * sqrt(1-[1]*[1]) + [2]",
    -0.5, 0.5);
  _circle_function_dn->SetParName(0, "radius");
  _circle_function_dn->SetParName(1, "sin(alpha)");
  _circle_function_dn->SetParName(2, "target");
}

//********************************************************************
bool InclinedTracks::ProcessEvent(const TEvent* event) {
//********************************************************************
  for (uint trackId = 0; trackId < event->GetTracks().size(); ++trackId) {

    TTrack* track = event->GetTracks()[trackId];
    if (!track)
      continue;

    if (_verbose > 1)
      std::cout << "Track id = " << trackId << std::endl;

    if (!sel::CrossingTrackSelection(track, _invert, _verbose))
      continue;

    std::vector<double> fit_v = sel::GetFitParams(track, _invert);
    std::vector<double> fit_xz = sel::GetFitParamsXZ(track, _invert);

    _angle_xy = fit_v[2];
    _angle_yz = fit_xz[2] * sel::v_drift_est;

    _ev = event->GetID();


    // diagonalise track
    std::vector<std::vector<THit*>> track_diag;
    for (auto col:track->GetCols(_invert)) {
      for (auto pad:col) {
        auto col_id = pad->GetCol(_invert);
        auto row_id = pad->GetRow(_invert);
        auto cons = col_id - row_id;
        if (col_id < 1 || col_id > 34 || row_id < 1 || row_id > 30)
          continue;
        if (cons >= 31 | cons <= -30)
          continue;
        // if (abs(cons%2) == 1)
        //   continue;

        // search if the diagonal is already considered
        std::vector<std::vector<THit*>>::iterator it;
        for (it = track_diag.begin(); it < track_diag.end(); ++it) {
          if (!(*it)[0]) {
            std::cerr << "error1" << std::endl;
            exit(1);
          }

          if ((*it)[0]->GetCol(_invert) -  (*it)[0]->GetRow(_invert) == cons) {
            (*it).push_back(pad);
            break;
          }
        } // loop over track_diag
        if (it == track_diag.end()) {
          std::vector<THit*> v;
          v.push_back(pad);
          track_diag.push_back(v);
        }
      } // over pads
    } // over cols diagonalise track

    // diag check
    // std::cout << "event" << std::endl;
    // for (auto col:track_diag) {
    //   std::cout << "column" << std::endl;
    //   for (auto pad:col) {
    //     std::cout << pad->GetCol(_invert) << "\t" << pad->GetRow(_invert) << "\t" << pad->GetCol(_invert) - pad->GetRow(_invert) << std::endl;
    //   }
    // }
    // std::cout << "event id " << event->GetID() << std::endl;
    // std::cout << "col mult "  <<  track_diag.size() << std::endl;

    // selection
    if (track_diag.size() < 25 )
      continue;
    _store_event = true;

    sort(track_diag.begin(), track_diag.end(), [&](std::vector<THit*> col1, std::vector<THit*> col2){
                    return GetInclinedX(col1[0]) < GetInclinedX(col2[0]);
                    });

    double track_pos[70];
    double cluster_mean[70];
    float charge_max[70];
    double a_peak_fit[70];

    // reset tree values
    for (auto colId = 0; colId < 70; ++colId) {
      _multiplicity[colId]  = -999;
      _charge[colId]        = -999;
      _residual[colId]      = -999;
      _clust_pos[colId]     = -999;
      _track_pos[colId]     = -999;
      _x[colId]             = -999;
      _x_av[colId]          = -999;
      _cluster_av[colId]    = -999;
      for (auto padId = 0; padId < 10; ++padId) {
        _dx[colId][padId]   = -999;
        _time[colId][padId]   = -999;
        _qfrac[colId][padId] = -999;
      }

      // cluster[colId]       = 0;
      // cluster_N[colId]     = 0;
      charge_max[colId]    = 0;
      track_pos[colId]     = -999.;
      cluster_mean[colId]  = 0.;
      a_peak_fit[colId]    = 0.;
    }

    // first loop over column
    std::vector<std::pair<float, float> > track_pos_vec;

    for (uint pairIt = 0; pairIt < track_diag.size(); pairIt += 2) {

      if (track_diag.size() - pairIt < 2)
        continue;

      std::pair<float, float> cluster[2];
      float cluster_e[2];
      for (uint pair = 0; pair < 2; ++pair) {
        auto colIt = pairIt + pair;
        auto col = track_diag[colIt];
        if (!col[0])
          continue;

        _multiplicity[colIt] = col.size();

        // centre of charge as fit starting point
        Float_t weight_mean = 0.;
        _charge[colIt] = 0;
        for (auto pad:col) {
          weight_mean += pad->GetQ() * GetInclinedY(pad);
          _charge[colIt] += pad->GetQ();
        }
        weight_mean /= _charge[colIt];

        track_pos[colIt] = weight_mean;

        // std::cout << "ini\t"  <<  track_pos[colIt] << std::endl;

        // fit the cluster with PRF
        auto chi2Function_cluster = [&](const Double_t *par) {
          //minimisation function computing the sum of squares of residuals
          // looping at the graph points
          double chi2 = 0;

          for (auto pad:col) {
            auto q      = pad->GetQ();
            if (!q)
              continue;

            double a = 1. * q / _charge[colIt];
            double center_pad_y = GetInclinedY(pad);

            // avoid using pads wich are far away from track
            // limit by PRF fitting range (PRF function robustness)
            // if (abs(center_pad_y - pos) > _fit_bound)
            //   continue;

            double part = (a - _PRF_function->Eval(par[0] - center_pad_y));
            if (_charge_uncertainty) {
              double c = 1.*_charge[colIt];
              double b = 1.*q;
              part *= c*c;
              part /= c*sqrt(b) + b*sqrt(c);
            }
            part *= part;

            chi2 += part;
          }
          return chi2;
        };

        if (_multiplicity[colIt] > 1) {
          ROOT::Math::Functor fcn_cluster(chi2Function_cluster,1);
          ROOT::Fit::Fitter  fitter_cluster;

          double pStart[1] = {track_pos[colIt]};
          fitter_cluster.SetFCN(fcn_cluster, pStart);
          fitter_cluster.Config().ParSettings(0).SetName("y");

          bool ok = fitter_cluster.FitFCN();
          (void)ok;
          const ROOT::Fit::FitResult & result_cluster = fitter_cluster.Result();

          _clust_pos[colIt] = result_cluster.GetParams()[0];
        } else {
          _clust_pos[colIt] = GetInclinedY(col[0]);
        }

        _x[colIt] = GetInclinedX(col[0]);
        cluster[pair] = std::make_pair(GetInclinedX(col[0]), _clust_pos[colIt]);
        if (_multiplicity[colIt] == 1)
          cluster_e[pair] = 0.01;
        else
          cluster_e[pair] = 0.008;

      } // loop over pair
      float av_x = 0.5*(cluster[0].first + cluster[1].first);
      float y1 = cluster[0].second;
      float y2 = cluster[1].second;
      float av_y = y1 * cluster_e[0] * cluster_e[0] + y2 * cluster_e[1] * cluster_e[1];
      av_y /= cluster_e[0] * cluster_e[0] + cluster_e[1] * cluster_e[1];
      track_pos_vec.push_back(std::make_pair(av_x, av_y));

    } // 1st loop over column

    // fit the track
    TGraphErrors* track_gr = new TGraphErrors();
    for (auto dot:track_pos_vec) {
      track_gr->SetPoint(track_gr->GetN(), dot.first, dot.second);
    }

    TF1* fit;
    TString func;
    if (!_do_arc_fit) {
      track_gr->Fit("pol1", "Q");
      fit = (TF1*)track_gr->GetFunction("pol1")->Clone();
    } else {
      Float_t q_up, q_down;
      q_up = q_down = 1.e9;
      _circle_function_dn->SetParameters(80., 0, 0.);
      track_gr->Fit("circle_dn", "Q");
      fit = track_gr->GetFunction("circle_dn");
      if (fit)
        q_down = fit->GetChisquare() / fit->GetNDF();

      _circle_function_up->SetParameters(80., 0, 0.);
      track_gr->Fit("circle_up", "Q");
      fit = track_gr->GetFunction("circle_up");
      if (fit)
        q_up = fit->GetChisquare() / fit->GetNDF();

      if (q_up > q_down)
        func = "circle_dn";
      else
        func = "circle_up";

      track_gr->Fit(func);
      if (!track_gr->GetFunction(func))
        continue;
      fit = (TF1*)track_gr->GetFunction(func)->Clone();
    }
    // end of track fit

    int colIt = 0;
    for (auto dot:track_pos_vec) {
      _x_av[colIt] = dot.first;
      double track_fit_y    = fit->Eval(dot.first);
      _cluster_av[colIt] = dot.second;
      _track_pos[colIt] = track_fit_y;
      _residual[colIt] = _cluster_av[colIt] - track_fit_y;
      _resol_col_hist[colIt]->Fill(_cluster_av[colIt] - track_fit_y);

      ++colIt;
    } // end of second loop

    for (colIt = 0; colIt < 70; ++colIt) {
      if (_clust_pos[colIt] == -999)
        continue;

      auto col = track_diag[colIt];
      int padId = 0;
      // std::cout << "column" << std::endl;
      for (auto pad:col) {
        auto inc_y = GetInclinedY(pad);
        auto inc_x = GetInclinedX(pad);
        auto q = pad->GetQ();
        auto time = pad->GetTime();

        double track_fit_y    = fit->Eval(inc_x);
        _dx[colIt][padId]      = inc_y - track_fit_y;
        _qfrac[colIt][padId]   = 1.*q / _charge[colIt];
        _time[colIt][padId]    = time;
        _PRF_histo->Fill( _dx[colIt][padId],
                          _qfrac[colIt][padId]);

        ++padId;
      }
    }
  } // loop over tracks

  if (_store_event) {
    if (!_batch)
      Draw();
    _tree->Fill();
    _passed_events.push_back(event->GetID());
  }

  return true;
}

bool InclinedTracks::Draw() {
  std::cout << "Draw event " << _ev << std::endl;
  TCanvas c("c_inc", "c_inc", 800, 600);
  TGraphErrors* gr = new TGraphErrors();
  TGraphErrors* gr_f = new TGraphErrors();
  TGraphErrors* gr_c = new TGraphErrors();
  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_cluster_av[colIt] == -999)
      continue;

    gr->SetPoint(gr->GetN(), 1e3*_x_av[colIt], 1e3*_cluster_av[colIt]);
    gr_f->SetPoint(gr_f->GetN(), 1e3*_x_av[colIt], 1e3*_track_pos[colIt]);
  }

  for (auto colIt = 0; colIt < 70; ++colIt) {
    if (_clust_pos[colIt] == -999)
      continue;

    gr_c->SetPoint(gr_c->GetN(), 1e3*_x[colIt], 1e3*_clust_pos[colIt]);
  }

  gr_c->SetTitle("Event " + TString().Itoa(_ev, 10));
  gr_c->SetMarkerStyle(kPlus);
  gr_c->Draw("ap");
  gr->Draw("same p");
  gr->GetYaxis()->SetTitle("Inclined Y, [mm]");
  gr->GetXaxis()->SetTitle("Inclined X, [mm]");

  gr_f->SetLineColor(kRed);
  gr_f->Draw("same l");
  c.Update();
  c.WaitPrimitive();
  return true;
}

Float_t InclinedTracks::GetInclinedX(const THit* h) {
  // 45 degree
  Float_t angle = 0.78539816;
  return geom::GetXpos(h->GetCol(_invert)) * TMath::Cos(angle) -
         geom::GetYpos(h->GetRow(_invert)) * TMath::Sin(angle);
}

Float_t InclinedTracks::GetInclinedY(const THit* h) {
  // 45 degree
  Float_t angle = 0.78539816;
  return   geom::GetXpos(h->GetCol(_invert)) * TMath::Sin(angle) +
           geom::GetYpos(h->GetRow(_invert)) * TMath::Cos(angle);
}























































