#include "DBSCANReconstruction.hxx"
#include <TCanvas.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>

DBSCANReconstruction::DBSCANReconstruction(): ReconstructionBase() {
  ;
}

bool DBSCANReconstruction::Initialize(int verbose) {
  std::cout << "Initialize DBSCAN Reconstruction.........";
  _verbose = verbose;
  std::cout << "done" << std::endl;
  return true;
}

double DBSCANReconstruction::MeasureDistance(Node a, Node b){
    double distance2 = pow(a.hit->GetRow() - b.hit->GetRow(), 2) +
                       pow(a.hit->GetCol()-b.hit->GetCol(),2) +
                       pow( (a.hit->GetTime()-b.hit->GetTime())/30 ,2); // distance squared.
    return pow(distance2,0.5);
}

std::vector<Node> DBSCANReconstruction::FindClusters(std::vector<Node> nodes){
  std::vector <Node> nodesToCheck;
  int clusterID = -1;
  int clustered = 0;
  while (clustered < (int) nodes.size()){
    if (!nodesToCheck.size()){
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
    while(nodesToCheck.size() > 0){
      int close_nodes = 0;
      Node origin = nodesToCheck[0];
      for (auto neighbor:nodes){
        if(origin.id == neighbor.id) continue;
        if(MeasureDistance(origin,neighbor) <= MIN_DIST){
          close_nodes++;
          if(neighbor.c < 0){
            bool found = std::any_of(nodesToCheck.begin(), nodesToCheck.end(),
                                [&](const Node& n) { return n.id == neighbor.id;}) ? true : false;
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

std::vector<Node> DBSCANReconstruction::FillNodes( TEvent* event){
  std::vector<Node> nodes;

  if (_verbose > 2)
    std::cout << "Found " << event->GetHits().size() << " hits" << std::endl;

  if (event->GetHits().size() > MAX_NODES_TOT)
    return nodes;

  // use in the future: instead of tmax, take all local maximums in the WF!!

  // int window = 20;
  // for(int i=0; i<geom::nPadx; i++){
  //   for(int j=0; j<geom::nPady; j++){
  //     for(int k=0; k<geom::Nsamples; k++){
  //       int Q = padAmpl[i][j][k];
  //       int maxNum = 0;
  //       bool maxFound = true;
  //       if(!Q) continue;
  //       for(int l=-window; l<=window; l++){
  //         if(k+l <0 or k+l >511) continue;
  //         if(padAmpl[i][j][k+l]>Q){
  //           maxFound = false;
  //           break;
  //         }
  //       }
  //       if(maxFound){
  //         maxNum++;
  //         Node node;
  //         node.x  = i;
  //         node.y  = j;
  //         node.t  = k;
  //         node.q  = padAmpl[i][j][k];
  //         node.id = nodes.size();
  //         nodes.push_back(node);
  //       }
  //     }
  //   }
  // }

  // fill the maximum amplitude and time
  for (auto hit:event->GetHits()) {
    if (hit->GetQ() > 0) {
      Node node;
      node.hit  = hit;
      node.id = nodes.size();
      nodes.push_back(node);
    }
  }

  return nodes;
}

// std::vector<int> DBSCANReconstruction::FillWFs(const TEvent* event, Node n) {
//   std::vector<int> wf_v;

//   for(int k=0; k<geom::Nsamples; k++)
//     wf_v.push_back(padAmpl[n.x][n.y][k]);

//   return wf_v;
// }

std::vector<DB_Cluster> DBSCANReconstruction::FindClustersLargerThan(std::vector<Node> nodes, int minNodes){
  std::vector<DB_Cluster> clusters;
  auto it_max = std::max_element(nodes.begin(), nodes.end(),
                       [](const Node& n1, const Node& n2) { return n1.hit->GetQ() < n2.hit->GetQ(); });
  int numClusters = (it_max == nodes.end()) ? -1 : (*it_max).c;
  int acceptedClusters = 0;
  for(int i=0; i<=numClusters; i++){
    DB_Cluster cluster;
    cluster.id = acceptedClusters;
    for (auto n:nodes) if (n.c == i){
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
    for(uint i = 0; i < clusters.size(); ++i)
      std::cout << "\t" << clusters[i].size;
    std::cout << std::endl;
  }
  return clusters;
}

std::vector<Node> DBSCANReconstruction::UpdateNodes(std::vector<DB_Cluster> clusters,
                                                    std::vector <Node> nodes){
    std::vector<Node> updated_nodes;
    for(auto n:nodes){
      n.c = -1;
      updated_nodes.push_back(n);
    }
    for(auto c:clusters){
        for(int i=0; i<c.size;i++){
            nodes[c.nodes[i]].c = c.id;
            updated_nodes.push_back(nodes[c.nodes[i]]);
        }
    }
    if (_verbose > 2)
      std::cout << "Updated " << updated_nodes.size() << " nodes" << std::endl;
    return updated_nodes;
}

void DBSCANReconstruction::DrawNodes(std::vector<Node> nodes){
  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.05);
  TH2F    *MM      = new TH2F("MM","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TH2F    *MMsel   = new TH2F("MMsel","",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");

  for (auto n:nodes){
    event3D->Fill(n.hit->GetTime(),n.hit->GetRow(),n.hit->GetCol(),n.hit->GetQ());
    MM->Fill(n.hit->GetCol(),n.hit->GetRow(),n.hit->GetQ());
    if(n.hit->GetQ()==0) MMsel->Fill(n.hit->GetCol(),n.hit->GetRow(),n.hit->GetQ());
  }

  TCanvas *canv = new TCanvas("canv", "canv", 800, 600, 800, 600);
  canv->Divide(3,1);
  canv->cd(1);
  MM->Draw("COLZ");
  canv->cd(2);
  MMsel->Draw("COLZ");

  canv->cd(3);
  event3D->Draw("x:y:z:c","","box2");
  TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetLimits(0,geom::nPadx);
  htemp->GetYaxis()->SetLimits(0,geom::nPady);
  htemp->GetZaxis()->SetLimits(0,500);
  htemp->SetTitle("");
  canv->Update();
  canv->WaitPrimitive();
  delete htemp;
  delete canv;

  delete MM;
  delete MMsel;
  delete event3D;
}

bool DBSCANReconstruction::FillOutput(TEvent* event,
                                      std::vector<Node> nodes,
                                      std::vector<DB_Cluster> clusters
                                      ){
  if(!nodes.size()) return false;
  std::vector<THit*> unusedHits;
  std::vector<int> usedHits;
  usedHits.resize(nodes.size(),0);
  // store selected hits
  for(uint trkID=0; trkID<clusters.size(); trkID++){
    if (!nodes.size())
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
    Node n = nodes[i];
    // std::vector<int> wf_v = FillWFs(event, n);
    if(!usedHits[i]) event->AddUnusedHit(n.hit);
  }

  return true;
}

bool DBSCANReconstruction::SelectEvent(TEvent* event) {
  auto raw_nodes = FillNodes(event);
  if (!raw_nodes.size())
    return false;
  std::vector<Node> nodes = FindClusters(raw_nodes);
  std::vector<DB_Cluster> clusters = FindClustersLargerThan(nodes,15);
  std::vector<Node> new_nodes = UpdateNodes(clusters,nodes);

  //if(nodes.size()) DrawNodes(new_nodes);
  if(nodes.size()) if(clusters.size() == 1) return FillOutput(event,
                                                              new_nodes,
                                                              clusters);
  return false;
}
