#include "TF1.h"
#include "TH2F.h"

#include "TrackSelection.hxx"

//******************************************************************************
bool TrackSel::CrossingTrackSelection(const TClusterPtrVec &track,
                                      const int  &max_mult,
                                      const float  &max_mean_mult,
                                      const bool &cut_gap,
                                      const float &max_phi,
                                      const float &max_theta,
                                      const std::vector<std::pair<int, int>>& _broken_pads,
                                      const bool &invert,
                                      const int  &verbose) {
//******************************************************************************
  Float_t m_mean;
  Int_t m_max;
  TrackSel::GetMultiplicity(track, m_mean, m_max);
  auto no_gap = TrackSel::GetNoGap(track, _broken_pads, invert);
  if (verbose > 1) {
    std::cout << "SELECTION " << std::endl;
    std::cout << "Max mult\t" << m_max << " < " << max_mult << std::endl;
    std::cout << "Mean mult\t" << m_mean << " < " << max_mean_mult << std::endl;
    std::cout << "No gap\t" << no_gap << std::endl;
    std::cout << "Linear Phi\t" << GetLinearPhi(track, invert) << std::endl;
    std::cout << "Linear theta\t" << GetLinearTheta(track, invert) << std::endl;
  }
  if (m_max > max_mult) return false;
  if (m_mean > max_mean_mult) return false;
  if (!no_gap && cut_gap) return false;

  if (max_phi > 0 && abs(GetLinearPhi(track, invert)) > max_phi) return false;
  if (max_theta > 0 && abs(GetLinearTheta(track, invert)) > max_theta) return false;

  return true;
}

//******************************************************************************
void TrackSel::GetMultiplicity(const TClusterPtrVec &track,
                               Float_t& m_mean,
                               Int_t& m_max
) {
//******************************************************************************
  m_max = -1;
  m_mean = 0.;
  Int_t n = 0;
  for (const auto& col:track) {
    auto mult = col->GetSize();
    if (mult < 1)
      continue;
    m_mean += num::cast<float>(mult);
    n += 1;
    if ((int)mult > m_max)
      m_max = num::cast<int>(mult);
  }
  m_mean /= num::cast<float>(n);
}

//******************************************************************************
bool TrackSel::GetNoGap(const TClusterPtrVec &track,
                        const std::vector<std::pair<int, int>>& _broken_pads,
                        const bool &invert
                   ) {
//******************************************************************************
  for (const auto & cluster:track) if (cluster->GetSize()) {
    std::vector<int> row;
    std::vector<int> col;
    for (const auto& pad:*cluster) if (pad) {
      if (pad->GetQ() < CHARGE_THR_FOR_GAP)
        continue;
      // broken pad fix --> always include the broken one
      // if the adjacent pad is there --> broken pad will not cause a gap
      for (const auto& broken : _broken_pads) {
        if (abs(pad->GetCol() - broken.first) < 2 && abs(pad->GetRow() - broken.second) < 2) {
          row.push_back(broken.second);
          col.push_back(broken.first);
        }
      }
      row.push_back(pad->GetRow());
      col.push_back(pad->GetCol());
    } // loop over pads

    sort(row.begin(), row.end());
    sort(col.begin(), col.end());

    // no gaps in rows
    if (!GetNoGapVector(row))
      return false;
    // no gaps in cols
    if (!GetNoGapVector(col))
      return false;
  } // loop over cluster

  return true;
}

bool TrackSel::GetNoGapVector(const std::vector<int>& v) {
  if (v.empty())
    return false;
  auto prev = v[0];
  for (const auto& r:v) {
    if (r - prev > 1)
      return false;
    prev = r;
  }
  return true;
}

//******************************************************************************
double TrackSel::GetLinearPhi(const TClusterPtrVec &track,
                              bool invert) {
//******************************************************************************
  std::vector <double> par = TrackSel::GetFitParams(track, invert);
  return par[2];
}

//******************************************************************************
double TrackSel::GetLinearTheta(const TClusterPtrVec &track,
                                bool invert) {
//******************************************************************************
  std::vector <double> par = TrackSel::GetFitParamsXZ(track, invert);
  return par[2] * TrackSel::v_drift_est;
}

//******************************************************************************
std::vector <double> TrackSel::GetFitParams(const TClusterPtrVec &track,
                                            bool invert) {
//******************************************************************************
  std::vector <double> params;
  params.reserve(3);
  for (auto i = 0; i < 3; ++i) params.push_back(-999.);

  TH2F* MM = new TH2F("MM", "MM",
                      geom::nPadx, 0, geom::nPadx,
                      geom::nPadx, 0, geom::nPadx
                      );

  for(auto & cluster:track)
    for (const auto& pad:*cluster)
      if(pad->GetQ()) {
        if (!invert)
          MM->Fill(pad->GetCol(),pad->GetRow(),pad->GetQ());
        else
          MM->Fill(pad->GetRow(),pad->GetCol(),pad->GetQ());
    }

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  if (fit){
    double quality = fit->GetChisquare() / fit->GetNDF();
    double b = fit->GetParameter(0);
    double k = fit->GetParameter(1);
    params[0] = quality;
    params[1] = b;
    params[2] = k;
  }

  delete MM;
  return params;
}


//******************************************************************************
std::vector <double> TrackSel::GetFitParamsXZ(const TClusterPtrVec &track,
                                              bool invert) {
//******************************************************************************
  std::vector <double> params;
  params.reserve(3);
  for (auto i = 0; i < 3; ++i) params.push_back(-999);

  TH2F* MM = new TH2F("MM", "MM",
                      geom::nPadx, 0, geom::nPadx,
                      geom::Nsamples, 0, geom::Nsamples
                      );
  for(auto & cluster:track) {
    auto q_lead = 0;
    auto x_lead = 0;
    auto y_lead = 0;
    auto z_lead = 0;
    for (const auto& pad:*cluster) if (pad->GetQ()) {
      if (pad->GetQ() > q_lead) {
        q_lead = pad->GetQ();
        x_lead = pad->GetCol();
        y_lead = pad->GetRow();
        z_lead = pad->GetTime();
      }
    }

    if (q_lead) {
      if (!invert)
        MM->Fill(x_lead, z_lead);
      else
        MM->Fill(y_lead, z_lead);
    }
  }

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  if (fit) {
    auto quality = fit->GetChisquare() / fit->GetNDF();
    double b = fit->GetParameter(0);
    double k = fit->GetParameter(1);
    params[0] = quality;
    params[1] = b;
    params[2] = k;
  }

  delete MM;
  return params;
}
