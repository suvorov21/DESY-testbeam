#include "Math/Functor.h"
#include "Fit/Fitter.h"

#include "TrackFitter.hxx"
#include "Geom.hxx"

//******************************************************************************
TrackFitterBase::TrackFitterBase() {
//******************************************************************************
  _circle_function = new TF1("circle",
    "-TMath::Sign(1, [0]) * sqrt(1./[0]*1./[0] - TMath::Power(x+0.198 - [1] * 1./[0], 2)) + 1./[0] * sqrt(1-[1]*[1]) + [2]",
    -0.5, 0.5);
  _circle_function->SetParName(0, "radius");
  _circle_function->SetParName(1, "sin(alpha)");
  _circle_function->SetParName(2, "target");
}

//******************************************************************************
std::pair<Double_t, Double_t> TrackFitCern::FitCluster(const THitPtrVec& col,
                                  const double& pos) {
//******************************************************************************
//  auto t_leading = -1;
  auto maxQ = -1;
  // WARNING pads should be already sorted
  auto sum = 0;
  for (auto & pad:col) {
    sum += pad->GetQ();
    if (pad->GetQ() > maxQ) {
      maxQ = pad->GetQ();
//      t_leading = pad->GetTime();
    }
  }

  auto chi2Function_cluster = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;

    for (auto & pad:col) {
      auto q      = pad->GetQ();
      auto colId  = pad->GetCol(_invert);
      if (!q)
        continue;

      double a = 1. * q / sum;
      double center_pad_y;
      center_pad_y = geom::GetYposPad(pad, _invert, _angle);

      // avoid using pads which are far away from track
      // limit by PRF fitting range (PRF function robustness)
      if (abs(center_pad_y - pos) > _fit_bound)
        continue;

      const TF1* PRF_tmp = _prf_function;

      if (_prf_size && _complicated_pattern_PRF) {
        auto prf_it = abs(colId % _prf_size);
        PRF_tmp = _prf_arr[prf_it];
      }

      if (_prf_size && _individual_column_PRF)
        PRF_tmp = _prf_arr[colId];

      if (!PRF_tmp) {
        throw std::logic_error("PRF is NULL");
      }

      double part = (a - PRF_tmp->Eval(center_pad_y - par[0]));

      if (_charge_uncertainty) {
        /** Poisson fluctuations only */
        double c = 1.*sum;
        double b = 1.*q;
        auto error_1 = c*sqrt(b) + b*sqrt(c);
        error_1 /= c*c;
        /** */

        /** Pedestal driven approach
        * sigma_q = 9; sigma_Q = 9*sqrt(n)
        */
        // auto error_2 = 1.*sum * 9 + 1.* sqrt(col.size()) * 9 * q;
        // error_2 /= 1.* sum * sum;
        /** */

        /** Profile uncertainty */
        // auto error_3 = 0.;
        // auto tol   = 999.;
        // for (auto i = 0; i < _PRF_graph->GetN(); ++i) {
        //   auto y = _PRF_graph->GetY()[i];
        //   if (abs(y - a) < tol) {
        //     tol   = abs(y - a);
        //     error_3 = _PRF_graph->GetEY()[i];
        //   }
        // }
        /** */

        part /= error_1;
      }
      part *= part;
      chi2 += part;

      /** Take time into account */
      // auto time = pad->GetTime();
      // if (t_leading <= 0 ||
      //     time - t_leading <= 0 ||
      //     time - t_leading > _PRF_time_func->GetMaximum() ||
      //     time - t_leading < _PRF_time_func->GetMinimum())
      //   continue;

      // auto x_time = abs(_PRF_time_func->GetX(1.*(time - t_leading)));
      // auto dt = _PRF_time_error->GetBinContent(_PRF_time_error->FindBin(time - t_leading));

      // double time_chi = abs(par[0] - center_pad_y) - x_time;
      // time_chi /= dt;
      // time_chi *= time_chi;

      // if (isinf(time_chi) || time_chi != time_chi)
      //   continue;

      // chi2 += time_chi;
      /** */
    }
    return chi2;
  };

  ROOT::Math::Functor fcn_cluster(chi2Function_cluster,1);
  ROOT::Fit::Fitter  fitter_cluster;

  double pStart[1] = {pos};
  fitter_cluster.SetFCN(fcn_cluster, pStart);
  fitter_cluster.Config().ParSettings(0).SetName("y");

  bool ok = fitter_cluster.FitFCN();
  (void)ok;
  const ROOT::Fit::FitResult & result_cluster = fitter_cluster.Result();

  return {result_cluster.GetParams()[0], result_cluster.GetErrors()[0]};
}

//******************************************************************************
std::shared_ptr<TF1> TrackFitCern::FitTrack(const TClusterPtrVec& clusters,
                                            const int& miss_id) {
//******************************************************************************
TGraphErrors track_gr;

  for (uint clusterId = 0; clusterId < clusters.size(); ++clusterId) {
    double x   = clusters[clusterId]->GetX();
    double y   = clusters[clusterId]->GetY();
    double y_e = clusters[clusterId]->GetYE();

    if (int(clusterId) == miss_id || isnan(x) || isnan(y))
      continue;

    if (_verbose > 3)
      std::cout << "Track point\t" << x << "\t" << y << std::endl;

    track_gr.SetPoint(track_gr.GetN(), x, y);
    track_gr.SetPointError(track_gr.GetN()-1, 0., y_e);
  } // loop over x

  if (_verbose > 2)
    std::cout << "Fit graph with " << track_gr.GetN() << " points" << std::endl;

  TF1* fit = nullptr;
  TString func;
  TString opt = "Q";
  if (_verbose > 2)
    opt = "";
  if (_shape == TrackShape::parabola) {
    track_gr.Fit("pol2", opt);
    fit = (TF1*)track_gr.GetFunction("pol2")->Clone();
  } else if (_shape == TrackShape::linear) {
    track_gr.Fit("pol1", opt);
    fit = (TF1*)track_gr.GetFunction("pol1")->Clone();
  } else if (_shape == TrackShape::arc) {
    // Float_t q_up, q_down;
    // q_up = q_down = 1.e9;
    // _circle_function_dn->SetParameters(80., 0, 0.);
    // track_gr.Fit("circle_dn", opt);
    // fit = track_gr.GetFunction("circle_dn");
    // if (fit)
    //   q_down = fit->GetChisquare() / fit->GetNDF();

    // _circle_function_up->SetParameters(80., 0, 0.);
    // track_gr.Fit("circle_up", opt);
    // fit = track_gr.GetFunction("circle_up");
    // if (fit)
    //   q_up = fit->GetChisquare() / fit->GetNDF();

    // if (q_up > q_down)
    //   func = "circle_dn";
    // else
    //   func = "circle_up";

    // track_gr.Fit(func, "Q");
    // if (!track_gr.GetFunction(func))
    //   return NULL;
    // fit = (TF1*)track_gr.GetFunction(func)->Clone();

    Double_t q_up, q_down;
    _circle_function->SetParameters(1./80., 0, 0.);
    track_gr.Fit("circle", opt);
    TF1* fit_up = (TF1*)track_gr.GetFunction("circle")->Clone();
    if (!fit_up)
      return nullptr;
    q_up = fit_up->GetChisquare() / fit_up->GetNDF();

    _circle_function->SetParameters(-1./80., 0, 0.);
    track_gr.Fit("circle", opt);
    TF1* fit_dn = (TF1*)track_gr.GetFunction("circle")->Clone();
    if (!fit_dn)
      return nullptr;
    q_down = fit_dn->GetChisquare() / fit_dn->GetNDF();

    if (!track_gr.GetFunction("circle"))
      return nullptr;

    if (q_up > q_down)
      fit = (TF1*)fit_dn->Clone();
    else
      fit = (TF1*)fit_up->Clone();

    delete fit_up;
    delete fit_dn;
  }

  if (!fit)
    return nullptr;

  return std::make_shared<TF1>(*fit);
}

//******************************************************************************
void TrackFitCern::SetPRFarr(TF1* f[], const int n) {
//******************************************************************************
  _prf_arr = new TF1*[n];
  for (auto i = 0; i < n; ++i) {
    _prf_arr[i] = f[i];
  }

  _prf_size = n;
}