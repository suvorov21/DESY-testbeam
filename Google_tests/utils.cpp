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

TEST(cast, normalCast) {
  int i = 5;
  auto f = num::cast<float>(i);
  EXPECT_EQ(f, i);
}

TEST(cast, errorThrowSign) {
  uint i = std::numeric_limits<int>::max() + 2;
  EXPECT_THROW(num::cast<int>(i), std::bad_cast);
}

TEST(cast, errorThrowPrecision) {
  int i = std::numeric_limits<int>::max();
  EXPECT_THROW(num::cast<float>(i), std::bad_cast);
}

TEST(cast, enoughPrecision) {
  int i = std::numeric_limits<int>::max();
  EXPECT_EQ(num::cast<double>(i), i);
}
