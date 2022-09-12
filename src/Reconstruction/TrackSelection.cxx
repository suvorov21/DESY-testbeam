#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"

#include "TrackSelection.hxx"

// 25 MHz --> 40 ns/bin   7 cm /us  -->   0.007 cm/ns ---> 0.28 cm / bin
// 50 Mhz --> ... --> 0.14 cm / bin
const float TrackSel::v_drift_est = 0.0028;

void TrackSel::Reset() {
    phi_ = 999;
    theta_ = 999;
}

//******************************************************************************
bool TrackSel::CrossingTrackSelection(const TClusterPtrVec &track) {
//******************************************************************************
    Float_t m_mean{10};
    Int_t m_max{10};
    GetMultiplicity(track, m_mean, m_max);
    auto no_gap = GetNoGap(track);

    FitXY(track);
    FitXZ(track);

    auto min_max_time = GetMinMaxTime(track);

    if (verbose_ > 1) {
        std::cout << "SELECTION " << std::endl;
        std::cout << "Max mult\t" << m_max << " < " << max_mult_ << std::endl;
        std::cout << "Mean mult\t" << m_mean << " < " << max_mean_mult_ << std::endl;
        std::cout << "No gap\t" << no_gap << std::endl;
        std::cout << "Linear Phi\t" << phi_ << std::endl;
        std::cout << "Linear theta\t" << theta_ << std::endl;
        std::cout << "Min max time\t" << min_max_time.first << "\t" << min_max_time.second << std::endl;
    }
    if (m_max > max_mult_) return false;
    if (m_mean > max_mean_mult_) return false;
    if (!no_gap && cut_gap_) return false;

    if (time_min_ > 0) {
        if (time_min_ > min_max_time.first) {
            return false;
        }
    }

    if (time_max_ > 0) {
        if (time_max_ < min_max_time.second) {
            return false;
        }
    }

    if (max_phi_ > 0 && abs(phi_) > max_phi_) return false;
    if (max_theta_ > 0 && abs(theta_) > max_theta_) return false;

    return true;
}

//******************************************************************************
void TrackSel::GetMultiplicity(const TClusterPtrVec &track,
                               Float_t &m_mean,
                               Int_t &m_max
) {
//******************************************************************************
    m_max = -1;
    m_mean = 0.;
    Int_t n = 0;
    for (const auto &col : track) {
        auto mult = col->GetSize();
        if (mult < 1)
            continue;
        m_mean += num::cast<float>(mult);
        n += 1;
        if ((int) mult > m_max)
            m_max = num::cast<int>(mult);
    }
    m_mean /= num::cast<float>(n);
}

//******************************************************************************
bool TrackSel::GetNoGap(const TClusterPtrVec &track) {
//******************************************************************************
    for (const auto &cluster : track)
        if (cluster->GetSize()) {
            std::vector<int> row;
            std::vector<int> col;
            for (const auto &pad : *cluster)
                if (pad) {
                    if (pad->GetQMax() < CHARGE_THR_FOR_GAP)
                        continue;
                    // broken pad fix --> always include the broken one
                    // if the adjacent pad is there --> broken pad will not cause a gap
                    for (const auto &broken : broken_pads_) {
                        if (pad->GetCard() == broken[0] && \
                            abs(pad->GetCol() - broken[1]) < 2 && \
                            abs(pad->GetRow() - broken[2]) < 2) {
                            row.push_back(broken[2]);
                            col.push_back(broken[1]);
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

bool TrackSel::GetNoGapVector(const std::vector<int> &v) {
    if (v.empty())
        return false;
    auto prev = v[0];
    for (const auto &r : v) {
        if (r - prev > 1)
            return false;
        prev = r;
    }
    return true;
}

//******************************************************************************
void TrackSel::FitXY(const TClusterPtrVec &track) {
//******************************************************************************
    TGraph gr;
    for (const auto& cluster : track) {
        gr.SetPoint(gr.GetN(),
                    Geom::GetXposPad((*cluster)[0], invert_),
                    Geom::GetYposPad((*cluster)[0], invert_)
                    );
    }
    gr.Fit(fit_->GetName(), "Q");
    auto fitFunc = gr.GetFunction(fit_->GetName());

    if (!fitFunc)
        return;

    double xStart = gr.GetX()[0];
    phi_ = TMath::Sin(TMath::ATan(fitFunc->Derivative(xStart)));
}

//******************************************************************************
void TrackSel::FitXZ(const TClusterPtrVec &track) {
//******************************************************************************
    TGraph gr;
    for (const auto& cluster : track) {
        gr.SetPoint(gr.GetN(),
                    Geom::GetXposPad((*cluster)[0], invert_),
                    (*cluster)[0]->GetTimeMax()
        );
    }
    gr.Fit(fit_->GetName(), "Q");
    auto fitFunc = gr.GetFunction(fit_->GetName());

    if (!fitFunc)
        return;

    double xStart = gr.GetX()[0];
    theta_ = TMath::Sin(TMath::ATan(fitFunc->Derivative(xStart) * TrackSel::v_drift_est));
}

//******************************************************************************
std::pair<int, int> TrackSel::GetMinMaxTime(const TClusterPtrVec &track) {
//******************************************************************************
    int min_time = 1000;
    int max_time = -1;
    for (const auto& cluster : track) {
        auto time = (*cluster)[0]->GetTimeMax();
        if (time > max_time) {
            max_time = time;
        }

        if (time < min_time) {
            min_time = time;
        }
    }
    return {min_time, max_time};
}