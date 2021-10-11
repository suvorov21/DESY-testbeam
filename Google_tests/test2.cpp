//
// Created by SERGEY SUVOROV on 10/10/2021.
//
#include "gtest/gtest.h"
#include "SpatialResolAna.hxx"
#include "AnalysisBase.hxx"

TEST(Initialisation, outputWrite) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->WriteOutput(), false);
}
