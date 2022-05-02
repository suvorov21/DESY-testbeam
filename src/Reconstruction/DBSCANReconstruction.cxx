#include "DBSCANReconstruction.hxx"
#include <TStyle.h>
#include <TF1.h>

DBSCANReconstruction::DBSCANReconstruction() : ReconstructionBase() {}

bool DBSCANReconstruction::Initialize(int verbose) {
    std::cout << "Initialize DBSCAN Reconstruction.........";
    _verbose = verbose;
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

    // FIXME put this in the selection
    if (module.size() > MAX_NODES_TOT)
        return nodes;

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
            if (pattern.second.size() > MIN_NODES_PER_PATT) {
                event->AddPattern(pattern.second);
            }
        }
    } // over modules

    MatchModules(event);

    return true;
}

void DBSCANReconstruction::MatchModules(const std::shared_ptr<TEvent> &event) {
    // FIXME
    // dummy algo
    // works for one module only
    for (const auto& pattern : event->GetAllPatterns()) {
        for (const auto& patternInMM : pattern.second) {
            event->AddTrack({patternInMM});
        }
    }

}
