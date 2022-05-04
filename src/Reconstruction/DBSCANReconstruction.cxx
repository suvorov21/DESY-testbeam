#include "DBSCANReconstruction.hxx"

#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"
#include "TVector3.h"

DBSCANReconstruction::DBSCANReconstruction() : ReconstructionBase() {}

bool DBSCANReconstruction::Initialize(int verbose) {
    std::cout << "Initialize DBSCAN Reconstruction.........";
    _verbose = verbose;
    trackFitFunction_ = std::make_unique<TF1>("TrackFit", "[0]+[1]*x+[2]*x*x", -3000., 0.);
    std::cout << "done" << std::endl;
    return true;
}

double DBSCANReconstruction::MeasureDistance(const Node &a, const Node &b) {
    double distance2 = pow(a.hit->GetRow() - b.hit->GetRow(), 2) +
        pow(a.hit->GetCol() - b.hit->GetCol(), 2) +
        pow((a.hit->GetTimeMax() - b.hit->GetTimeMax()) / 30, 2); // distance squared.
    return pow(distance2, 0.5);
}

std::vector<Node> DBSCANReconstruction::FindClusters(std::vector<Node> &nodes) {
    std::vector<Node> nodesToCheck;
    int clusterID = -1;
    int clustered = 0;
    while (clustered < num::cast<int>(nodes.size())) {
        if (nodesToCheck.empty()) {
            clusterID++;
            auto it = std::find_if(nodes.begin(), nodes.end(),
                                   [](const Node &n) { return n.c < 0; });
            if (it != nodes.end()) {
                nodesToCheck.push_back(*(it));
                nodes[(*it).id].c = clusterID;
                clustered++;
            }
        }
        std::vector<Node> tmpNodes;
        while (!nodesToCheck.empty()) {
            int close_nodes = 0;
            Node origin = nodesToCheck[0];
            for (auto neighbor : nodes) {
                if (origin.id == neighbor.id) continue;
                if (MeasureDistance(origin, neighbor) <= MIN_DIST) {
                    close_nodes++;
                    if (neighbor.c < 0) {
                        bool found =
                            std::any_of(nodesToCheck.begin(), nodesToCheck.end(),
                                        [&](const Node &n) { return n.id == neighbor.id; });
                        if (!found) tmpNodes.push_back(neighbor);
                    }
                }
            }
            if (close_nodes >= MIN_NODES) {
                nodesToCheck.insert(nodesToCheck.end(), tmpNodes.begin(), tmpNodes.end());
                if (nodes[origin.id].c < 0) {
                    clustered++;
                    nodes[origin.id].c = clusterID;
                }
            }
            tmpNodes.clear();
            nodesToCheck.erase(nodesToCheck.begin());
        }
    }
    if (_verbose > 2)
        std::cout << "Found " << nodes.size() << " nodes in " << clusterID << " clusters" << std::endl;
    return nodes;
}

std::vector<Node> DBSCANReconstruction::FillNodes(const THitPtrVec &module) {
    std::vector<Node> nodes;

    if (_verbose > 2)
        std::cout << "Found " << module.size() << " hits" << std::endl;

    // fill the maximum amplitude and time
    for (const auto &hit : module) {
        if (hit->GetQMax() > 0) {
            Node node;
            node.hit = hit;
            node.id = num::cast<int>(nodes.size());
            nodes.push_back(node);
        }
    }

    return nodes;
}

bool DBSCANReconstruction::ReconstructEvent(const std::shared_ptr<TEvent> &event) {
    auto hitsMap = event->GetAllHits();
    // pattern recognition per modules
    for (auto const& module : hitsMap) {
        auto raw_nodes = FillNodes(module.second);
        if (raw_nodes.empty())
            continue;
        std::vector<Node> nodes = FindClusters(raw_nodes);

        // loop over reconstructed clusters
        std::unordered_map<int, THitPtrVec> patterns{};
        for (const auto& node : nodes) {
            if (node.c >= 0)
                patterns[node.c].emplace_back(node.hit);
        }

        for (auto const& pattern : patterns) {
            if (pattern.second.size() > MIN_HITS_PER_PATT) {
                event->AddPattern(pattern.second);
            }
        }
    } // over modules

    MatchModules(event);

    return true;
}

void DBSCANReconstruction::MatchModules(const std::shared_ptr<TEvent> &event) {

    std::vector<std::pair<short, std::ptrdiff_t>> groveStone;

    std::vector<TPatternVec> fTracks{};

    for (short mmStart = 0; mmStart < Geom::nModules; ++mmStart) {
        // for each MM determine the neighbors
        auto neighbors = getNeighboursMM(mmStart);
        // loop over them
        for (const auto mmSecond : neighbors) {
            // mmSecond is always following mmStart
            // if not --> the possible merge has been already done
            if (mmSecond < mmStart)
                continue;

            // for each pair of trajectories in the adjacent MMs
            auto patternsInFirstModule = event->GetPatternsInModule(mmStart);
            for (auto trajFirstIt  = patternsInFirstModule.begin();
                      trajFirstIt != patternsInFirstModule.end();
                      ++trajFirstIt) {
                auto trajFirst = *trajFirstIt;
                // loop over trajes in the adjacent module
                for (auto& trajSecond : event->GetPatternsInModule(mmSecond)) {
                    if (fitTogether(trajFirst, trajSecond)) {
                        // match found
                        groveStone.emplace_back(mmStart, std::distance(patternsInFirstModule.begin(), trajFirstIt));
                        // add hits from traj 1 to traj2
                        for (const auto& hit : trajFirst)
                            trajSecond.emplace_back(hit);
                    } // if fit together
                } // over 2nd Traj
            } // over 1st Traj
        } // over second MM
    } // over first MM

    // Fill the summary info
    for (short mmStart = 0; mmStart < Geom::nModules; ++mmStart) {
        auto patternsInFirstModule = event->GetPatternsInModule(mmStart);
        for (auto trajFirstIt = patternsInFirstModule.begin();
             trajFirstIt != patternsInFirstModule.end();
             ++trajFirstIt) {

            // if the pattern was not merged later
            if (std::find_if(groveStone.cbegin(),
                             groveStone.cend(),
                             [&](const auto &p1) {
                               return p1.first == mmStart && \
                          p1.second == std::distance(patternsInFirstModule.begin(),
                                                     trajFirstIt);
                             }) != groveStone.cend())
                 continue;
            // split patterns back into modules
            std::unordered_map<short, TPattern> splitToModules;
            for (const auto& hit : *trajFirstIt) {
                splitToModules[hit->GetCard()].emplace_back(hit);
            }

            TTrack track;
            for (const auto& pattern : splitToModules) {
                track.emplace_back(pattern.second);
            }
            event->AddTrack(track);
        } // over trajs
    } // over modules

    // dummy algo
    // works for one module only
//    for (const auto& pattern : event->GetAllPatterns()) {
//        for (const auto& patternInMM : pattern.second) {
//            event->AddTrack({patternInMM});
//        }
//    }

}


std::vector<short> DBSCANReconstruction::getNeighboursMM(const ushort i) {
    std::vector<short> v;
    // see Geom.hxx for layout
    if (i != 0 && i != 4) {
        v.emplace_back(i-1);
        // diagonal
        v.emplace_back(i + (i < 4 ? 3 : -5));
    }
    if (i != 3 && i != 7) {
        v.emplace_back(i+1);
        // diagonal
        v.emplace_back(i + (i < 4 ? 5 : -3));
    }
    if (i > 3) {
        v.emplace_back(i-4);
    }
    if (i < 4) {
        v.emplace_back(i+4);
    }
    return v;
}

bool DBSCANReconstruction::fitTogether(const TPattern &traj1, const TPattern &traj2) {
    TGraph traj1Graph;
    TGraph traj2Graph;
    TGraph joinedGraph;

    // FIXME consider doing it one time per traj. Right now it's done for every pair of trajs
    for (const auto &hit : traj1) {
        if (!hit)
            continue;
        if (hit->GetQMax() < chargeThresholdFit_)
            continue;
        traj1Graph.SetPoint(traj1Graph.GetN(),
                            Geom::GetXposPad(hit),
                            Geom::GetYposPad(hit));
        joinedGraph.SetPoint(joinedGraph.GetN(),
                             Geom::GetXposPad(hit),
                             Geom::GetYposPad(hit));
    } // loop over hits

    for (const auto &hit : traj2) {
        if (!hit)
            continue;
        if (hit->GetQMax() < chargeThresholdFit_)
            continue;

        traj2Graph.SetPoint(traj2Graph.GetN(),
                            Geom::GetXposPad(hit),
                            Geom::GetYposPad(hit));
        joinedGraph.SetPoint(joinedGraph.GetN(),
                             Geom::GetXposPad(hit),
                             Geom::GetYposPad(hit));
    } // loop over hits

    // fit trajs separately and together
    traj1Graph.Fit(trackFitFunction_->GetName(), "Q");
    traj2Graph.Fit(trackFitFunction_->GetName(), "Q");
    joinedGraph.Fit(trackFitFunction_->GetName(), "Q");

    auto fit1 = traj1Graph.GetFunction(trackFitFunction_->GetName());
    auto fit2 = traj2Graph.GetFunction(trackFitFunction_->GetName());
    auto fitJ = joinedGraph.GetFunction(trackFitFunction_->GetName());
    // fit fails
    if (!fit1 || !fit2 || !fitJ) {
        if (_verbose > 0)
            std::cout << "One of fit fails\n";
        return false;
    }

    auto Q1 = fit1->GetChisquare() / fit1->GetNDF();
    auto Q2 = fit2->GetChisquare() / fit2->GetNDF();
    auto QJ = fitJ->GetChisquare() / fitJ->GetNDF();

    if (_verbose > 0)
        std::cout << "Fit quality: " << Q1 << "\t" << Q2 << "\t" << QJ << std::endl;

    // if any of individual fits is "bad"
    if (Q1 > fitQthreshold_ || Q2 > fitQthreshold_ || QJ > fitQthreshold_)
        return false;

    // if join fit is good comparing to separate fits
    if (QJ < sqrt(Q1 * Q2) * 1.2 || QJ < std::max(Q1, Q2))
        return true;

    return false;
}
