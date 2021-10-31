#include "TF1.h"
#include "TH2F.h"

#include "Selection.hxx"

//******************************************************************************
bool sel::CrossingTrackSelection( const std::vector<std::unique_ptr<TCluster>> &track,
                                  const int  &max_mult,
                                  const bool &cut_gap,
                                  const float &max_phi,
                                  const float &max_theta,
                                  const bool &invert,
                                  const int  &verbose) {
//******************************************************************************

  auto m_max = sel::GetMaxMultiplicity(track);
  auto no_gap = sel::GetNoGap(track, invert);
  if (verbose > 1) {
    std::cout << "SELECTION " << std::endl;
    std::cout << "Max mult\t" << m_max << " < " << max_mult << std::endl;
    std::cout << "No gap\t" << no_gap << std::endl;
    std::cout << "Linear Phi\t" << GetLinearPhi(track, invert) << std::endl;
    std::cout << "Linear theta\t" << GetLinearTheta(track, invert) << std::endl;
  }
  if (m_max > max_mult) return false;
  if (!no_gap && cut_gap) return false;

  if (max_phi > 0 && abs(GetLinearPhi(track, invert)) > max_phi) return false;
  if (max_theta > 0 && abs(GetLinearTheta(track, invert)) > max_theta) return false;

  return true;
}

//******************************************************************************
int sel::GetMaxMultiplicity(const std::vector<std::unique_ptr<TCluster>> &track) {
//******************************************************************************
  auto it_max = std::max_element(track.begin(), track.end(),
                       [](const std::unique_ptr<TCluster> & cl1,
                                const std::unique_ptr<TCluster> & cl2) {
                                  return cl1->GetSize() < cl2->GetSize();
                        });
  if (it_max != track.end())
    return (int)(*it_max)->GetSize();
  return 999;
}

//******************************************************************************
bool sel::GetNoGap(const std::vector<std::unique_ptr<TCluster>> &track,
                   const bool &invert) {
//******************************************************************************
  for(auto & cluster:track) if(cluster->GetSize()){
    std::vector<int> row;
    std::vector<int> col;
    for (auto pad:*cluster) if (pad) {
      row.push_back(pad->GetRow(invert));
      col.push_back(pad->GetCol(invert));
    } // loop over pads

    sort(row.begin(), row.end());
    sort(col.begin(), col.end());

    // no gaps in rows
    if (row.empty())
      return false;
    auto prev = row[0];
    for (auto r:row) {
      if (r - prev > 1)
        return false;
      prev = r;
    }
    // no gaps in cols
    prev = col[0];
    for (auto c:col) {
      if (c - prev > 1)
        return false;
      prev = c;
    }
  } // loop over cluster

  return true;
}

//******************************************************************************
double sel::GetLinearPhi(const std::vector<std::unique_ptr<TCluster>> &track,
                        bool invert) {
//******************************************************************************
  std::vector <double> par = sel::GetFitParams(track, invert);
  return par[2];
}

//******************************************************************************
double sel::GetLinearTheta(const std::vector<std::unique_ptr<TCluster>> &track,
                          bool invert) {
//******************************************************************************
  std::vector <double> par = sel::GetFitParamsXZ(track, invert);
  return par[2] * sel::v_drift_est;
}

//******************************************************************************
std::vector <double> sel::GetFitParams(const std::vector<std::unique_ptr<TCluster>> &track,
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
    for (auto pad:*cluster)
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
std::vector <double> sel::GetFitParamsXZ(const std::vector<std::unique_ptr<TCluster>> &track,
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
    for (auto pad:*cluster) if (pad->GetQ()) {
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
