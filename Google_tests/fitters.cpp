//
// Created by SERGEY SUVOROV on 28/11/2021.
//

#include "gtest/gtest.h"

#include "FakeMaker.h"
#include "SpatialResolAna.hxx" // for PRF initialisation
#include "TrackFitter.hxx"

TEST(FitterTest, exepFitterBase) {
  auto shape = TrackFitterBase::arc;
  auto fitter = new TrackFitterBase(shape, false, 1, 0);
  EXPECT_THROW(fitter->FitCluster(), std::logic_error);
  EXPECT_THROW(fitter->FitTrack(), std::logic_error);
}

TEST(FitterTest, exepEmptyPRF) {
  auto fitter = std::make_unique<TrackFitCern>();
  std::vector<std::shared_ptr<THit>> col;
  col.emplace_back(std::make_shared<THit>());
  col[0]->SetQ(10);
  EXPECT_THROW(fitter->FitCluster(col), std::logic_error);
}

TEST(FitterTest, ClusterFit) {
  auto func = SpatialResolAna::InitializePRF("prf");

  auto fitter = std::make_unique<TrackFitCern>(TrackFitterBase::arc, false, 0, 1, func, nullptr, 0.027, true, nullptr, nullptr, 0);
  auto center = fake::GetCenteredCluster();
  auto result = fitter->FitCluster(center, geom::GetYposPad(center[0]) + 0.003);
  EXPECT_NEAR(result, geom::GetYposPad(center[0]), 0.001);

  center = fake::GetSideCLuster();
  result = fitter->FitCluster(center, geom::GetYposPad(center[0]));
  EXPECT_NEAR(result, geom::GetYposPad(center[0]) + 0.005, 0.001);
}
