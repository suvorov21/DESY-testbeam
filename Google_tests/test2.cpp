//
// Created by SERGEY SUVOROV on 10/10/2021.
//
#include "gtest/gtest.h"
//#include "SpatialResolAna.hxx"
#include "AnalysisBase.hxx"

#include "TRawEvent.hxx"

#include "Geom.hxx"

TEST(Initialisation, outputWrite) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->WriteOutput(), false);
}

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