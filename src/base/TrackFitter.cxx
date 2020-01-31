#include "Math/Functor.h"
#include "Fit/Fitter.h"

#include "TrackFitter.hxx"
#include "Geom.hxx"

//********************************************************************
TrackFitter::TrackFitter(FitterType type, TF1* func,
  float fit_bound, float uncertainty, Int_t it, Int_t verbose, bool invert,
  bool charge_uncertainty, bool do_arc_fit) {
  //********************************************************************
  _type         = type;
  _PRF_function = func;

  _fit_bound    = fit_bound;
  _uncertainty  = uncertainty;

  _iteration  = it;
  _verbose    = verbose;
  _invert     = invert;

  _charge_uncertainty = charge_uncertainty;
  _do_arc_fit = do_arc_fit;

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
Double_t TrackFitter::FitCluster(const std::vector<THit*>& col,
    const int cluster, double pos, pads_t& pos_in_pad,
     TH1F* uncertainty) {
  //********************************************************************

  switch (_type) {
    case CERN_like:     return GetClusterPosCERN(col, cluster, pos);
    case ILC_like:      return GetClusterPosILC(col, pos);
    case Separate_pads: return GetClusterPosInPad(col, cluster,
                                                  pos, pos_in_pad,
                                                  uncertainty); break;
    std::cout << "ERROR. TrackFitter::FitCluster() ";
    std::cout << "uknown type " << _type << std::endl;
    exit(1);
  }

  return 0.;

}

//********************************************************************
TF1* TrackFitter::FitTrack(const double* track_pos, const int* mult,
  const TTrack* track, const double pos, const pads_t pos_in_pad,
  const int miss_id) {
  //********************************************************************
  switch (_type) {
    case CERN_like:     return GetTrackFitCERN(track_pos, mult, miss_id);
    case ILC_like:      return GetTrackFitILC(track, pos, miss_id);
    case Separate_pads: return GetTrackFitSeparatePad(pos_in_pad, miss_id);

    std::cout << "ERROR. TrackFitter::FitTrack() ";
    std::cout << "uknown type " << _type << std::endl;
    exit(1);
  }

  return NULL;
}

//********************************************************************
double TrackFitter::GetClusterPosInPad(const std::vector<THit*>& col,
    const int cluster, const double pos,
    pads_t& pos_in_pad,  TH1F* uncertainty) {
  //********************************************************************
  if (!_iteration) {
    pos_in_pad[col[0]->GetCol(_invert)].push_back(std::make_pair(0.1, std::make_pair(pos, default_error)));
    return pos;
  }

  double sum1 = 0;
  double sum2 = 0;

  double track_pos;
  double track_pos_err;

  for (auto pad:col) {
    auto q      = pad->GetQ();
    auto it_y   = pad->GetRow(_invert);
    auto it_x   = pad->GetCol(_invert);
    if (!q)
      continue;

    if (_verbose > 5)
      std::cout << "pad " << "\t" << 1.*q/cluster << std::endl;

    double center_pad_y = geom::GetYpos(it_y, _invert);
    if (1.*q/cluster < _PRF_function->Eval(_fit_bound))
      continue;

    if (_verbose > 5)
      std::cout << q << "\t" << cluster << std::endl;

    if (1.*q/cluster > _PRF_function->GetParameter(0)) {
      track_pos = center_pad_y;
      track_pos_err = 0.003;
    } else {
      double pos_bias = _PRF_function->GetX(1.*q/cluster, 0., _fit_bound);

      if (pos > center_pad_y)
        track_pos = center_pad_y + pos_bias;
      else
        track_pos = center_pad_y - pos_bias;

      track_pos_err = sigma_pedestal / cluster;
      track_pos_err /= abs(_PRF_function->Derivative(pos_bias));

      if (_verbose > 5) {
        std::cout << "errors\t" << sigma_pedestal / cluster << "\t";
        std::cout << _PRF_function->Derivative(pos_bias) << "\t";
        std::cout << track_pos_err << std::endl;
      }
    }

    if (1.*q/cluster > 0.5 && track_pos_err > 0.003)
      track_pos_err = 0.003;

    if (track_pos_err > 0.08)
      track_pos_err = 0.08;

    if (uncertainty)
      track_pos_err = uncertainty->GetBinContent(uncertainty->FindBin(1.*q/cluster));

    if (_verbose > 5) {
      std::cout << "pad pos \t" << track_pos << "\t";
      std::cout << track_pos_err << std::endl;
    }

    pos_in_pad[it_x].push_back(std::make_pair(1.*q/cluster,
      std::make_pair(track_pos, track_pos_err)));

    sum1 += track_pos * TMath::Power(track_pos_err, -2);
    sum2 += TMath::Power(track_pos_err, -2);
  } // loop over pads in cluster

  return sum1 / sum2;
}


//********************************************************************
double TrackFitter::GetClusterPosCERN(const std::vector<THit*>& col,
    const int cluster, const double pos) {
  //********************************************************************

  if (!_iteration)
    return pos;

  auto chi2Function_cluster = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;

      double a = 1. * q / cluster;
      double center_pad_y = geom::GetYpos(it_y, _invert);

      // avoid using pads wich are far away from track
      // limit by PRF fitting range (PRF function robustness)
      if (abs(center_pad_y - pos) > _fit_bound)
        continue;

      double part = (a - _PRF_function->Eval(par[0] - center_pad_y));
      if (_charge_uncertainty) {
        double c = 1.*cluster;
        double b = 1.*q;
        part *= c*c;
        part /= c*sqrt(b) + b*sqrt(c);
      }
      part *= part;

      chi2 += part;
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

  return result_cluster.GetParams()[0];
}

// BUG PRF limitation is not implemented in the ILC fitting procedure
// At the moment consider usage of the Gaus-Lorentzian w/o limitation
//********************************************************************
double TrackFitter::GetClusterPosILC(const std::vector<THit*>& col,
  const double pos) {
  //********************************************************************

  auto chi2Function_cluster = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;

    double a_nom = 0.;
    double a_den = 0.;
    double a_tot = 0.;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;
      double center_pad_y = geom::GetYpos(it_y, _invert);

      a_nom += _PRF_function->Eval(center_pad_y - par[0]);
      a_den += TMath::Power(_PRF_function->Eval(center_pad_y - par[0]), 2) / q;
    }
    a_tot = a_nom / a_den;

    for (auto pad:col) {
      auto q      = pad->GetQ();
      auto it_y   = pad->GetRow(_invert);
      if (!q)
        continue;
      double center_pad_y = geom::GetYpos(it_y, _invert);

      double part = (q - a_tot*_PRF_function->Eval(center_pad_y - par[0]));
      part *= part;
      part /= q;

      chi2 += part;
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

  return result_cluster.GetParams()[0];
}

//********************************************************************
TF1* TrackFitter::GetTrackFitCERN(const double* track_pos,
  const int* mult, const int miss_id) {
  //********************************************************************
  TGraphErrors* track_gr = new TGraphErrors();

  for (auto it_x = 0; it_x < geom::GetMaxColumn(_invert); ++it_x) {
    double x = geom::GetXpos(it_x, _invert);
    if (track_pos[it_x] == -999.)
      continue;

    if (it_x == miss_id)
      continue;

    track_gr->SetPoint(track_gr->GetN(), x, track_pos[it_x]);
    double error;
    if (mult[it_x] == 1)
      error = one_pad_error;
    else {
      if (_iteration == 0)
        error = default_error;
      else
        error = _uncertainty;
    }
    track_gr->SetPointError(track_gr->GetN()-1, 0., error);
  } // loop over x

  if (_verbose > 2)
    std::cout << "Fit graph with " << track_gr->GetN() << " points" << std::endl;

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

    track_gr->Fit(func, "Q");
    if (!track_gr->GetFunction(func))
      return NULL;
    fit = (TF1*)track_gr->GetFunction(func)->Clone();

  }

  delete track_gr;

  if (!fit)
    return NULL;

  return fit;
}

//********************************************************************
TF1* TrackFitter::GetTrackFitSeparatePad(const pads_t pos_in_pad,
  const int miss_id) {
  //********************************************************************
  TGraphErrors* track_gr = new TGraphErrors();

  for (auto it_x = 0; it_x < geom::GetMaxColumn(_invert); ++it_x) {
    double x = geom::GetXpos(it_x, _invert);
    if (!pos_in_pad[it_x].size())
      continue;

    if (it_x == miss_id)
      continue;

    for (auto it = pos_in_pad[it_x].begin(); it < pos_in_pad[it_x].end(); ++it) {
      track_gr->SetPoint(track_gr->GetN(), x, (*it).second.first);
      track_gr->SetPointError(track_gr->GetN()-1, 0., (*it).second.second);
    }
  } // loop over x

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

    track_gr->Fit(func, "Q");
    if (!track_gr->GetFunction(func))
      return NULL;
    fit = (TF1*)track_gr->GetFunction(func)->Clone();

  }

  delete track_gr;

  if (!fit)
    return NULL;

  return fit;
}

//********************************************************************
TF1* TrackFitter::GetTrackFitILC(const TTrack* track, const double pos,
  const int miss_id) {
  //********************************************************************
  auto chi2Function = [&](const Double_t *par) {
    //minimisation function computing the sum of squares of residuals
    // looping at the graph points
    double chi2 = 0;
    float y_pos = 0;

    if (par[3] > 0)
      _circle_function_up->SetParameters(par[0], par[1], par[2]);
    else
      _circle_function_dn->SetParameters(par[0], par[1], par[2]);


    for (auto col:track->GetCols(_invert)) {
      auto it_x = col[0]->GetCol(_invert);
      auto x    = geom::GetXpos(it_x, _invert);

      // TODO this is redefinition of the SpaceResolution::MissColumn()
      if (it_x == 0 || it_x == geom::GetMaxColumn(_invert)-1)
        continue;

      if (it_x == miss_id)
        continue;

      double a_nom = 0.;
      double a_den = 0.;
      double a_tot = 0.;

      if (par[3] > 0) {
        y_pos = _circle_function_up->Eval(x);
      } else {
        y_pos = _circle_function_dn->Eval(x);
      }

      for (auto pad:col) {
        auto q      = pad->GetQ();
        auto it_y   = pad->GetRow(_invert);
        if (!q)
          continue;

        double center_pad_y = geom::GetYpos(it_y, _invert);

        a_nom += _PRF_function->Eval(center_pad_y - y_pos);
        a_den += TMath::Power(_PRF_function->Eval(center_pad_y - y_pos), 2) / q;
      }
      a_tot = a_nom / a_den;

      for (auto pad:col) {
        auto q      = pad->GetQ();
        auto it_y   = pad->GetRow(_invert);
        if (!q)
          continue;
        double center_pad_y = geom::GetYpos(it_y, _invert);

        double part = (q - a_tot*_PRF_function->Eval(center_pad_y - y_pos));
        part *= part;
        part /= q;

        chi2 += part;
      }
    }
    return chi2;
  };

  float quality_dn, quality_up;
  quality_dn = quality_up = 1.e9;

  ROOT::Math::Functor fcn(chi2Function,4);
  ROOT::Fit::Fitter  fitter_dn;

  // fitting arc down
  double pStart[4] = {80., 0., pos, -1};
  fitter_dn.SetFCN(fcn, pStart, 4, true);
  fitter_dn.Config().ParSettings(0).SetName("R");
  fitter_dn.Config().ParSettings(1).SetName("sin(#alpha)");
  fitter_dn.Config().ParSettings(2).SetName("Target");
  fitter_dn.Config().ParSettings(3).SetName("Arc_dir (up/down)");

  fitter_dn.Config().ParSettings(0).SetLimits(0., 1.e5);
  fitter_dn.Config().ParSettings(3).Fix();

  bool ok = fitter_dn.FitFCN();
  (void)ok;
  const ROOT::Fit::FitResult & result_dn = fitter_dn.Result();

  _circle_function_dn->SetParameters( result_dn.GetParams()[0],
                                      result_dn.GetParams()[1],
                                      result_dn.GetParams()[2]);
  if (_verbose > 2)
    result_dn.Print(std::cout);

  quality_dn = result_dn.Chi2() / (track->GetHits().size() - 3);

  // fitting arc up
  ROOT::Math::Functor fcn_up(chi2Function,4);
  ROOT::Fit::Fitter  fitter_up;
  double pStart_up[4] = {80., 0., pos, 1};
  fitter_up.SetFCN(fcn_up, pStart_up, 4, true);
  fitter_up.Config().ParSettings(0).SetName("R");
  fitter_up.Config().ParSettings(1).SetName("sin(#alpha)");
  fitter_up.Config().ParSettings(2).SetName("Target");
  fitter_up.Config().ParSettings(3).SetName("Arc_dir (up/down)");
  fitter_up.Config().ParSettings(0).SetLimits(0., 1.e5);
  fitter_up.Config().ParSettings(3).Fix();

  ok = fitter_up.FitFCN();

  const ROOT::Fit::FitResult & result_up = fitter_up.Result();

  _circle_function_up->SetParameters( result_up.GetParams()[0],
                                      result_up.GetParams()[1],
                                      result_up.GetParams()[2]);
  if (_verbose > 2)
    result_up.Print(std::cout);

  quality_up = result_up.Chi2() / (track->GetHits().size() - 3);

  TF1* fit;
  TString func;
  float quality;
  if (quality_up > quality_dn) {
    fit = _circle_function_dn;
    quality = quality_dn;
  } else {
    fit = _circle_function_up;
    quality = quality_up;
  }

  if (quality > 100.)
    return NULL;

  //_chi2->Fill(quality);

  return fit;
}
