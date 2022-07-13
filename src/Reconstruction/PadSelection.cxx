//
// Created by SERGEY SUVOROV on 09/04/2022.
//

#include "PadSelection.hxx"
#include "Geom.hxx"

//******************************************************************************
THitPtrVec PadSelection::GetRobustPadsInCluster(THitPtrVec col,
                                                const broken_pads_t &broken_pads) {
//******************************************************************************
    std::vector<std::shared_ptr<THit>> result;
    // sort in charge decreasing order
    sort(col.begin(), col.end(), [](const std::shared_ptr<THit> &hit1,
                                    const std::shared_ptr<THit> &hit2) {
      return hit1->GetQMax() > hit2->GetQMax();
    });

    // leading pad
    auto moduleMM = col[0]->GetCard();
    auto col_id = col[0]->GetCol();
    auto row_id = col[0]->GetRow();
    // excluded from analysis the whole cluster if leading pad is near the broken pad
    for (const auto &broken : broken_pads) {
        if (moduleMM == broken[0] && \
            abs(col_id - broken[1]) < 2 && \
            abs(row_id - broken[2]) < 2)
            return result;
    }

    for (const auto &pad : col) {
        auto q = pad->GetQMax();
        if (!q)
            continue;

        /** cross-talk candidate */
        // if (i > 0 &&
        //     pad->GetTime() - col[0]->GetTime() < 4 &&
        //     1.*q / col[0]->GetQMax() < 0.08)
        //   continue;
        /** */

        // not more then 3 pads
        // if (i > 1)
        //   continue;

        // // WF with negative dt
        // if (pad->GetTime() - col[0]->GetTime() < -1)
        //   continue;

        // // avoid "suspicious" WF with small time difference in the 3rd pad
        // if (i > 1 && pad->GetTime() - col[0]->GetTime() < 5)
        //   continue;

        result.push_back(pad);

        // auto it_y   = pad->GetRow(_invert);
        // auto center_pad_y = geom::GetYpos(it_y, _invert);
    }

    return result;
}

//******************************************************************************
TClusterPtrVec PadSelection::GetRobustClusters(TClusterPtrVec &tr) {
//******************************************************************************
    TClusterPtrVec result;
    // sort clusters in increasing order
    sort(tr.begin(), tr.end(), [](TClusterPtr &cl1,
                                  TClusterPtr &cl) {
      return cl1->GetCharge() < cl->GetCharge();
    });

    // truncation cut
    /* NO TRUNCATION */
    auto frac = 1.00;
    auto i_max = num::cast<int>(round(frac * num::cast<double>(tr.size())));
    result.reserve(i_max);
    for (auto i = 0; i < i_max; ++i) {
        result.push_back(std::move(tr[i]));
    }
    /* */

    /* truncate with prominence */
    // sort along the track
    // sort(tr.begin(), tr.end(), [](TCluster* cl1,
    //                               TCluster* cl){
    //                                 return  cl1->GetX() < cl->GetX();});
    // compute the prominence
    // auto prom_cut = 0.6;
    // for (uint i = 1; i < tr.size() - 1; ++i) {
    //   auto prom = 2.*tr[i]->GetCharge() / (tr[i-1]->GetCharge() + tr[i+1]->GetCharge());
    //   if (prom > prom_cut)
    //     result.push_back(tr[i]);
    // }
    /* */

    // BUG truncation with neighbours is not working with clusters
    // trancation + neibours
    // auto frac = 0.95;
    // std::vector<Int_t> bad_pads;
    // Int_t i_max = round(frac * tr.size());
    // for (uint i = i_max; i < tr.size(); ++i) {
    //   bad_pads.push_back(tr[i][0]->GetCol(_invert));
    // }
    // for (auto i = 0; i < i_max; ++i) {
    //   auto it_x = tr[i][0]->GetCol(_invert);
    //   if (find(bad_pads.begin(), bad_pads.end(), it_x+1) == bad_pads.end() ||
    //       find(bad_pads.begin(), bad_pads.end(), it_x-1) == bad_pads.end())
    //     result.push_back(tr[i]);
    // }

    // cut on the total charge in the cluster
    // auto q_cut = 2000;
    // for (auto col:tr) {
    //   auto total_q = accumulate(col.begin(), col.end(), 0,
    //                       [](const int& x, const THit* hit)
    //                       {return x + hit->GetQMax();}
    //                       );
    //   if (total_q < q_cut)
    //     result.push_back(col);
    // }

    // sort by X for return
    sort(result.begin(), result.end(),
         [&](TClusterPtr &cl1, TClusterPtr &cl2) {
           return cl1->GetX() < cl2->GetX();
         });
    return result;
}

