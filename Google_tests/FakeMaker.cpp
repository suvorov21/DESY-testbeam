//
// Created by SERGEY SUVOROV on 28/11/2021.
//

#include "FakeMaker.h"
#include "Geom.hxx"
#include <array>

std::shared_ptr<TEvent> fake::GetEvent() {
  auto event = std::make_shared<TEvent>();
  event->SetID(2);
  auto hit = std::make_shared<THit>(10, 10, 200, 20);
  event->AddUsedHit(hit);
  return event;
}

THitPtrVec fake::GetTrack() {
  THitPtrVec track;
  for (auto col = 0; col < Geom::nPadx; ++col) {
    for (auto row = 15; row < 18; ++row) {
      auto hit = std::make_shared<THit>(col, row, 100, 200);
      track.push_back(hit);
    }
  }
  return track;
}

TClusterPtrVec fake::GetClusterdTrack() {
  TClusterPtrVec track;
  for (auto col = 0; col < Geom::nPadx; ++col) {
    THitPtrVec hits = GetCenteredCluster();
    TClusterPtr cluster = std::make_shared<TCluster>();
    for (const auto& hit:hits)
      cluster->AddHit(hit);

    track.emplace_back(cluster);
  }
  return track;
}

THitPtrVec fake::GetCenteredCluster() {
  THitPtrVec cluster;
  std::array<THitPtr, 3> hits;
  hits[0] = std::make_shared<THit>(10, 10, 0, 0, 0, 100, 800);
  hits[1] = std::make_shared<THit>(10, 11, 0, 0, 0, 150, 100);
  hits[2] = std::make_shared<THit>(10, 9, 0, 0, 0, 150, 100);
  std::copy(hits.begin(), hits.end(), std::back_inserter(cluster));
  return cluster;
}

THitPtrVec fake::GetSideCLuster() {
  THitPtrVec cluster;
  std::array<THitPtr, 3> hits;
  hits[0] = std::make_shared<THit>(10, 10, 0, 0, 0, 120, 600);
  hits[1] = std::make_shared<THit>(10, 11, 0, 0, 0, 120, 600);
  hits[2] = std::make_shared<THit>(10, 9, 0, 0, 0, 200, 700);
  std::copy(hits.begin(), hits.end(), std::back_inserter(cluster));
  return cluster;
}
