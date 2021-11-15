//
// Created by SERGEY SUVOROV on 10/10/2021.
//
#include "gtest/gtest.h"

#include "TrackFitter.hxx"
#include "TRawEvent.hxx"
#include "Geom.hxx"

TEST(GeomThrowTest, OferflowX) {
  EXPECT_THROW(geom::GetXpos(geom::GetNColumn() + 1),
              std::logic_error
               );
}

TEST(GeomThrowTest, OferflowY) {
  EXPECT_THROW(geom::GetYpos(geom::GetNColumn() + 1),
               std::logic_error
               );
}

TEST(CtorTests, hitCtor) {
  auto hit = std::make_unique<THit>();
  hit->SetADC(1, 10);
  hit->SetADC(1000, 10);

  ASSERT_EQ(hit->GetADC(1), 10);
  ASSERT_EQ(hit->GetADC(1000), 0);
}

TEST(CtorTests, eventCtor) {
  auto event = std::make_unique<TRawEvent>();
  ASSERT_EQ(event->GetID(), 0);
}

TEST(CtorTests, clusterCtor) {
  auto cluster = std::make_unique<TCluster>();
  ASSERT_EQ(cluster->GetCharge(), 0);
}

TEST(fitterTest, exepFitterBase) {
  auto shape = TrackFitterBase::arc;
  auto fitter = new TrackFitterBase(shape, false, 1, 0);
  EXPECT_THROW(fitter->FitCluster(), std::logic_error);
  EXPECT_THROW(fitter->FitTrack(), std::logic_error);
}

TEST(fitterTest, exepEmptyPRF) {
  auto fitter = std::make_unique<TrackFitCern>();
  std::vector<std::shared_ptr<THit>> col;
  col.emplace_back(std::make_shared<THit>());
  col[0]->SetQ(10);
  EXPECT_THROW(fitter->FitCluster(col, 0, 0), std::logic_error);
}