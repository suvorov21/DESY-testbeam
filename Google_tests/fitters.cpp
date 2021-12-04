//
// Created by SERGEY SUVOROV on 28/11/2021.
//

#include "gtest/gtest.h"

#include "FakeMaker.h"
#include "SpatialResolAna.hxx" // for PRF initialisation
#include "TrackFitter.hxx"

class Fitter: public ::testing::Test {
public:
  Fitter() {
    auto func = SpatialResolAna::InitializePRF("prf");
    fitter = std::make_unique<TrackFitCern>();
    fitter->SetPRF(func);
    fitter->SetFitBound(0.027);
    fitter->SetChargeUncertainty(true);
  }
  std::unique_ptr<TrackFitCern> fitter;
};

TEST(FitterTest, exepEmptyPRF) {
  auto fitter = std::make_unique<TrackFitCern>();
  std::vector<std::shared_ptr<THit>> col;
  col.emplace_back(std::make_shared<THit>());
  col[0]->SetQ(10);
  EXPECT_THROW(fitter->FitCluster(col, 0.), std::logic_error);
}

TEST_F(Fitter, ClusterFit) {
  auto center = fake::GetCenteredCluster();
  auto result = fitter->FitCluster(center, geom::GetYposPad(center[0]) + 0.003);
  EXPECT_NEAR(result.first, geom::GetYposPad(center[0]), 0.001);

  center = fake::GetSideCLuster();
  result = fitter->FitCluster(center, geom::GetYposPad(center[0]));
  EXPECT_NEAR(result.first, geom::GetYposPad(center[0]) + geom::dy/2, 0.001);
}

TEST_F(Fitter, TrackFit) {
  auto track = fake::GetClusterdTrack();
  for (auto& cluster : track) {
    auto result = fitter->FitCluster(cluster->GetHits(), geom::GetYposPad((*cluster)[0]));
    cluster->SetY(result.first);
  }
  auto fit = fitter->FitTrack(track, -1);
  EXPECT_NE(fit, nullptr);
  fitter->SetTrackShape(TrackFitterBase::TrackShape::arc);
  fit = fitter->FitTrack(track, -1);
  EXPECT_NE(fit, nullptr);
}
