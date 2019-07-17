#include "DBSCANSelection.hxx"
#include <TCanvas.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>

DBSCANSelection::DBSCANSelection(): SelectionBase() {
  ;
}

bool DBSCANSelection::Initialize() {
  std::cout << "Initialize crossing selection............";

  std::cout << "done" << std::endl;
  return true;
}

double DBSCANSelection::MeasureDistance(Node a, Node b){
    double distance2 = pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow( (a.t-b.t)/1 ,2); // distance squared.
    return pow(distance2,0.5);
}

std::vector<Node> DBSCANSelection::FindClusters(std::vector<Node> nodes){
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

std::vector<Node> DBSCANSelection::FillNodes(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples]){
  std::vector<Node> nodes;
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

std::vector<Cluster> DBSCANSelection::FindClustersLargerThan(std::vector<Node> nodes, int minNodes){
  std::vector<Cluster> clusters;
  int numClusters = 0;
  for (auto n:nodes) if(n.c > numClusters) numClusters = n.c;
  int acceptedClusters = 0;
  for(int i=0; i<numClusters; i++){
    Cluster cluster;
    cluster.id = acceptedClusters;
    for (auto n:nodes) if (n.c == i){
      cluster.nodes[cluster.size] = n.id;
      cluster.size++;
    }
    if(cluster.size >= minNodes){
      acceptedClusters++;
      clusters.push_back(cluster);
    }
  }
  // std::cout << "# of clusters: " << numClusters << std::endl; 
  // std::cout << "# of accepted clusters: " << acceptedClusters << std::endl; 
  return clusters;
}

std::vector<Node> DBSCANSelection::UpdateNodes(std::vector<Cluster> clusters, std::vector <Node> nodes){
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

void DBSCANSelection::DrawNodes(std::vector<Node> nodes){
  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.05);
  TH2F    *MM      = new TH2F("MM","MM",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");

  for (auto n:nodes){ 
    event3D->Fill(n.t,n.y,n.x,n.c);
    MM->Fill(n.x,n.y,n.q);
  }

  TCanvas *canv = new TCanvas("canv", "canv", 800, 600, 800, 600);
  canv->Divide(2,1);
  canv->cd(1);

  MM->Draw("COLZ"); 

  canv->cd(2);
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
  delete event3D;
}

bool DBSCANSelection::CheckQuality(std::vector<Node> nodes){
  TH2F    *MM      = new TH2F("MM","MM",geom::nPadx,0,geom::nPadx,geom::nPady,0,geom::nPady);
  TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");

  for (auto n:nodes){ 
    event3D->Fill(n.t,n.y,n.x,n.c);
    MM->Fill(n.x,n.y,n.q);
  }

  MM->Fit("pol1", "Q");
  TF1* fit = MM->GetFunction("pol1");

  double quality = 1.0e10;
  if (fit){
    quality = fit->GetChisquare() / fit->GetNDF();
    double k = fit->GetParameter(1);
    double b = fit->GetParameter(0);
    if(quality<1.0e6){
      delete MM;
      delete event3D;
      return true;
    }
  }

  delete MM;
  delete event3D;

  return false;
}

bool DBSCANSelection::FillOutput(std::vector<Node> nodes, int numTracks, Event& event){
  
  for(int num=0; num<numTracks; num++){  
    std::vector<std::vector<Int_t>> temp_evt = GetEmptyEvent();
    for (int itx = 0; itx < geom::nPadx; ++itx)
      for (uint ity = 0; ity < geom::nPady; ++ity){      
        temp_evt[itx][ity] = 0;
        // std::cout << itx << "," << ity << ", " << temp_evt[itx][ity] << std::endl; 
      }
    for(auto n:nodes){
      temp_evt[n.x][n.y] = n.q;
    }
    event.trackNum = numTracks;
    event.twoD.push_back(temp_evt);
  }

  // std::cout << "trackNum:  " << event.trackNum << std::endl;
  // std::cout << "eventsNum: " << event.twoD.size() << std::endl; 

  return true;
}

bool DBSCANSelection::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                  Event& event) {

  std::vector<Node> nodes = FindClusters(FillNodes(padAmpl));
  std::vector<Cluster> clusters = FindClustersLargerThan(nodes,15);
  std::vector<Node> new_nodes = UpdateNodes(clusters,nodes);
  // if(nodes.size()) DrawNodes(new_nodes);
  //if(nodes.size()>50) return CheckQuality(new_nodes);

  //if(nodes.size()) if(clusters.size() == 1) if(CheckQuality(new_nodes)) DrawNodes(new_nodes);
  if(nodes.size()) /*if(clusters.size() == 1)*/ if(CheckQuality(new_nodes)) return FillOutput(new_nodes,clusters.size(),event);
  return false;
}
