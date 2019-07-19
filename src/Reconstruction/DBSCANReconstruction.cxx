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

bool DBSCANReconstruction::Initialize() {
  std::cout << "Initialize DBSCAN Reconstruction............";

  std::cout << "done" << std::endl;
  return true;
}

double DBSCANReconstruction::MeasureDistance(Node a, Node b){
    double distance2 = pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow( (a.t-b.t)/30 ,2); // distance squared.
    return pow(distance2,0.5);
}

std::vector<Node> DBSCANReconstruction::FindClusters(std::vector<Node> nodes){
  std::vector <Node> nodesToCheck;
  int clusterID = -1;
  int clustered = 0;
  while (clustered < (int) nodes.size()){
    if (!nodesToCheck.size()){
      clusterID++;
      for(auto frst:nodes) if(frst.c < 0) {
        nodesToCheck.push_back(frst);
        nodes[frst.id].c = clusterID;
        clustered++;
        break;
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
            bool found = false;                         
            for (auto nTc:nodesToCheck){
              if(nTc.id == neighbor.id) found = true;
            }
            if(!found) tmpNodes.push_back(neighbor); 
          }
        }
      }
      if(close_nodes >= MIN_NODES){
        for(auto tmp:tmpNodes) nodesToCheck.push_back(tmp);
        if(nodes[origin.id].c<0){
          clustered++;
          nodes[origin.id].c = clusterID;
        }
      }
      tmpNodes.clear();
      nodesToCheck.erase(nodesToCheck.begin());
    }
  }
  return nodes;
}

std::vector<Node> DBSCANReconstruction::FillNodes(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples]){
  std::vector<Node> nodes;

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

  for(int i=0; i<geom::nPadx; i++){
    for(int j=0; j<geom::nPady; j++){
      int Q = 0;
      int t = 0;
      for(int k=0; k<geom::Nsamples; k++){
        if(!padAmpl[i][j][k]) continue;
        if(padAmpl[i][j][k]>Q) {Q = padAmpl[i][j][k]; t=k;}
      }
      if(Q){
        Node node;
        node.x  = i;
        node.y  = j;
        node.t  = t;
        node.q  = Q;
        node.id = nodes.size();
        nodes.push_back(node);
      }
    }
  }
  return nodes;
}

std::vector<Cluster> DBSCANReconstruction::FindClustersLargerThan(std::vector<Node> nodes, int minNodes){
  std::vector<Cluster> clusters;
  int numClusters = 0;
  int maxNodes = 150;
  for (auto n:nodes) if(n.c > numClusters) numClusters = n.c;
  int acceptedClusters = 0;
  for(int i=0; i<=numClusters; i++){
    Cluster cluster;
    cluster.id = acceptedClusters;
    for (auto n:nodes) if (n.c == i){
      cluster.nodes[cluster.size] = n.id;
      cluster.size++;
    }
    if(cluster.size >= minNodes && cluster.size <= maxNodes){
      acceptedClusters++;
      clusters.push_back(cluster);
    }
  }
  return clusters;
}

std::vector<Node> DBSCANReconstruction::UpdateNodes(std::vector<Cluster> clusters, std::vector <Node> nodes){
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
    event3D->Fill(n.t,n.y,n.x,n.c);
    MM->Fill(n.x,n.y,n.q);
    if(n.c==0) MMsel->Fill(n.x,n.y,n.q);
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

bool DBSCANReconstruction::FillOutput(std::vector<Node> nodes, std::vector<Cluster> clusters, Event& event){
  int numTracks = clusters.size();
  int trkCNT = 0;
  event.ResizeHits(nodes.size());
  for(int trkID=0; trkID<numTracks; trkID++){ 
    event.ResizeTracks(numTracks);
    event.tracks[trkCNT].ResizeCols();
    event.tracks[trkCNT].ResizeRows();
    event.tracks[trkCNT].ResizeHits(clusters[trkID].size);
    int sel_id = 0;
    int tot_id = 0;
    for(auto n:nodes){
      Hit hit;
      hit.c = n.x;
      hit.r = n.y;
      hit.t = n.t;
      hit.q = n.q;
      hit.id = tot_id;
      tot_id++;
      event.SetHit(hit.id,hit);
      if(n.c == trkCNT){
        hit.id = sel_id;
        sel_id++;
        event.tracks[trkCNT].SetHit(hit.id,hit);
        event.tracks[trkCNT].PushBackCol(hit.c,hit.id);
        event.tracks[trkCNT].PushBackRow(hit.r,hit.id);
      }
    }
    trkCNT++;
  } 
  return true;
}

bool DBSCANReconstruction::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                  Event& event) {

  std::vector<Node> nodes = FindClusters(FillNodes(padAmpl));
  std::vector<Cluster> clusters = FindClustersLargerThan(nodes,15);
  std::vector<Node> new_nodes = UpdateNodes(clusters,nodes);

  //if(nodes.size()) DrawNodes(new_nodes);
  if(nodes.size()) if(clusters.size() == 1) return FillOutput(new_nodes,clusters,event);
  return false;
}
