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
