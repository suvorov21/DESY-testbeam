#include "DBSCANReconstruction.hxx"
#include <TStyle.h>
#include <TF1.h>

DBSCANReconstruction::DBSCANReconstruction(): ReconstructionBase() {}

bool DBSCANReconstruction::Initialize(int verbose) {
  std::cout << "Initialize DBSCAN Reconstruction.........";
  _verbose = verbose;
  std::cout << "done" << std::endl;
  return true;
}

double DBSCANReconstruction::MeasureDistance(const Node& a, const Node& b){
    double distance2 = pow(a.hit->GetRow() - b.hit->GetRow(), 2) +
                       pow(a.hit->GetCol()-b.hit->GetCol(),2) +
                       pow( (a.hit->GetTime()-b.hit->GetTime())/30 ,2); // distance squared.
    return pow(distance2,0.5);
}

std::vector<Node> DBSCANReconstruction::FindClusters(std::vector<Node>& nodes){
  std::vector <Node> nodesToCheck;
  int clusterID = -1;
  int clustered = 0;
  while (clustered < (int) nodes.size()){
    if (nodesToCheck.empty()){
      clusterID++;
      auto it = std::find_if(nodes.begin(), nodes.end(),
                             [](const Node& n) {return n.c < 0;});
      if (it != nodes.end()) {
        nodesToCheck.push_back(*(it));
        nodes[(*it).id].c = clusterID;
        clustered++;
      }
    }
    std::vector <Node> tmpNodes;
    while(!nodesToCheck.empty()){
      int close_nodes = 0;
      Node origin = nodesToCheck[0];
      for (auto neighbor:nodes){
        if(origin.id == neighbor.id) continue;
        if(MeasureDistance(origin,neighbor) <= MIN_DIST){
          close_nodes++;
          if(neighbor.c < 0){
            bool found =
                std::any_of(nodesToCheck.begin(), nodesToCheck.end(),
                            [&](const Node &n) { return n.id == neighbor.id; });
            if(!found) tmpNodes.push_back(neighbor);
          }
        }
      }
      if(close_nodes >= MIN_NODES){
        nodesToCheck.insert(nodesToCheck.end(), tmpNodes.begin(), tmpNodes.end());
        if(nodes[origin.id].c<0){
          clustered++;
          nodes[origin.id].c = clusterID;
        }
      }
      tmpNodes.clear();
      nodesToCheck.erase(nodesToCheck.begin());
    }
  }
  if (_verbose > 2)
    std::cout << "Found " << nodes.size() << " nodes" << std::endl;
  return nodes;
}

std::vector<Node> DBSCANReconstruction::FillNodes(const std::shared_ptr<TEvent>& event){
  std::vector<Node> nodes;

  if (_verbose > 2)
    std::cout << "Found " << event->GetHits().size() << " hits" << std::endl;

  if (event->GetHits().size() > MAX_NODES_TOT)
    return nodes;

  // fill the maximum amplitude and time
  for (const auto& hit:event->GetHits()) {
    if (hit->GetQ() > 0) {
      Node node;
      node.hit  = hit;
      node.id = (int)nodes.size();
      nodes.push_back(node);
    }
  }

  return nodes;
}

std::vector<DB_Cluster> DBSCANReconstruction::FindClustersLargerThan(const std::vector<Node>& nodes, int minNodes){
  std::vector<DB_Cluster> clusters;
  auto it_max = std::max_element(nodes.begin(), nodes.end(),
                       [](const Node& n1, const Node& n2) { return n1.hit->GetQ() < n2.hit->GetQ(); });
  int numClusters = (it_max == nodes.end()) ? -1 : (*it_max).c;
  int acceptedClusters = 0;
  for(int i=0; i<=numClusters; i++){
    DB_Cluster cluster;
    cluster.id = acceptedClusters;
    for (const auto& n:nodes) if (n.c == i){
      // TODO make it a vector
      cluster.nodes[cluster.size] = n.id;
      cluster.size++;
    }
    if(cluster.size >= minNodes && cluster.size <= MAX_NODES){
      acceptedClusters++;
      clusters.push_back(cluster);
    }
  }
  if (_verbose > 2) {
    std::cout << "Found larger " << clusters.size() << " nodes";
    for(auto & cluster : clusters)
      std::cout << "\t" << cluster.size;
    std::cout << std::endl;
  }
  return clusters;
}

std::vector<Node> DBSCANReconstruction::UpdateNodes(const std::vector<DB_Cluster>& clusters,
                                                    std::vector <Node>& nodes){
    std::vector<Node> updated_nodes;
    for(auto n:nodes){
      n.c = -1;
      updated_nodes.push_back(n);
    }
    for(const auto & c:clusters){
        for(int i=0; i<c.size;i++){
            nodes[c.nodes[i]].c = c.id;
            updated_nodes.push_back(nodes[c.nodes[i]]);
        }
    }
    if (_verbose > 2)
      std::cout << "Updated " << updated_nodes.size() << " nodes" << std::endl;
    return updated_nodes;
}

bool DBSCANReconstruction::FillOutput(const std::shared_ptr<TEvent>& event,
                                      const std::vector<Node>& nodes,
                                      const std::vector<DB_Cluster>& clusters
                                      ){
  if(nodes.empty()) return false;
  std::vector<std::shared_ptr<THit>> unusedHits;
  std::vector<int> usedHits;
  usedHits.resize(nodes.size(),0);
  // store selected hits
  for(uint trkID=0; trkID<clusters.size(); trkID++){
    if (nodes.empty())
      continue;
    for (uint i = 0; i<nodes.size(); i++){
      Node n = nodes[i];
      if(n.c == (int)trkID){
        usedHits[i] = 1;
        event->AddUsedHit(n.hit);
      }
    }
  }

  // stored unselected hits
  for (uint i = 0; i<nodes.size(); i++){
    // std::vector<int> wf_v = FillWFs(event, n);
    if(!usedHits[i]) event->AddUnusedHit(nodes[i].hit);
  }

  return true;
}

bool DBSCANReconstruction::SelectEvent(const std::shared_ptr<TEvent>& event) {
  auto raw_nodes = FillNodes(event);
  if (raw_nodes.empty())
    return false;
  std::vector<Node> nodes = FindClusters(raw_nodes);
  std::vector<DB_Cluster> clusters = FindClustersLargerThan(nodes,15);
  std::vector<Node> new_nodes = UpdateNodes(clusters,nodes);

  //if(nodes.size()) DrawNodes(new_nodes);
  if(!nodes.empty()) if(clusters.size() == 1) return FillOutput(event,
                                                                new_nodes,
                                                                clusters);
  return false;
}
