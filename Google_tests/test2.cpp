//
// Created by SERGEY SUVOROV on 10/10/2021.
//
#include "gtest/gtest.h"
#include "SpatialResolAna.hxx"
#include "AnalysisBase.hxx"

#include "Geom.hxx"

TEST(Initialisation, outputWrite) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->WriteOutput(), false);
}

TEST(GeomDeathTest, OferflowX) {
  EXPECT_THROW(geom::GetXpos(geom::GetNColumn() + 1),
              std::logic_error
               );
}

TEST(GeomDeathTest, OferflowY) {
  EXPECT_THROW(geom::GetYpos(geom::GetNColumn() + 1),
               std::logic_error
               );
}