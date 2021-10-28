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
  auto hit = new THit();
  hit->SetADC(1, 10);
  hit->SetADC(1000, 10);

  ASSERT_EQ(hit->GetADC(1), 10);
  ASSERT_EQ(hit->GetADC(1000), 0);
}

TEST(CtorTests, eventCtor) {
  auto event = new TRawEvent();
  ASSERT_EQ(event->GetID(), 0);
}

TEST(CtorTests, clusterCtor) {
  auto cluster = new TCluster();
  ASSERT_EQ(cluster->GetCharge(), 0);
  delete cluster;
}

TEST(fitterTest, exepFitterBase) {
  auto shape = TrackFitterBase::arc;
  auto fitter = new TrackFitterBase(shape, false, 1, 0);
  EXPECT_THROW(fitter->FitCluster(), std::logic_error);
  EXPECT_THROW(fitter->FitTrack(), std::logic_error);
}

TEST(fitterTest, exepEmptyPRF) {
  auto fitter = new TrackFitCern();
  std::vector<THit*> col;
  col.emplace_back(new THit());
  col[0]->SetQ(10);
  EXPECT_THROW(fitter->FitCluster(col, 0, 0), std::logic_error);
}