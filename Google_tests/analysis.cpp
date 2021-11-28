//
// Created by SERGEY SUVOROV on 18/09/2021.
//
#include "gtest/gtest.h"
#include "AnalysisBase.hxx"
#include "DummyMaker.h"

TEST(InitialisationTest, ParamFileRead) {
  auto ana = std::make_unique<AnalysisBase>();
  EXPECT_EQ(ana->ReadParamFile(), true);
}

TEST(InitialisationTest, ParamFileNotRead) {
  auto ana = std::make_unique<AnalysisBase>();
  ana->setParamFile("blabla.ini");
  EXPECT_EQ(ana->ReadParamFile(), false);
}

TEST(InitialisationTest, CLIread) {
  auto argc = 2;
  char argv1[] = "main";
  char argv2[] = "test";
  char argv3[] = "-v4";

  char* argv_test1[2];
  argv_test1[0] =  argv1;
  argv_test1[1] =  argv2;

  char* argv_test2[2];
  argv_test2[0] =  argv1;
  argv_test2[1] =  argv3;

  auto ana = std::make_unique<AnalysisBase>();
  EXPECT_THROW(ana->ReadCLI(argc, argv_test1), std::logic_error);
  EXPECT_EQ(ana->ReadCLI(argc, argv_test2), true);
}

TEST(InitialisationTest, outputWriteErr) {
  auto ana = std::make_unique<AnalysisBase>();
  EXPECT_EQ(ana->WriteOutput(), false);
}

TEST(InitialisationTest, inputException) {
  auto ana = new AnalysisBase();
  ana->setParamFile("blabla.ini");
  EXPECT_DEATH(ana->Initialize(), "Parameter file is not read");
  ana->setParamFile("");
  EXPECT_DEATH(ana->Initialize(), "No input file specified");
  ana->setInputFile("input.root");
  EXPECT_DEATH(ana->Initialize(), "No output file specified");
  ana->setOutputFile("output.root");
  EXPECT_DEATH(ana->Initialize(), "Unknown tree name");
}

TEST(InitialisationTest, baseProcessEvent) {
  auto ana = std::make_unique<AnalysisBase>();
  auto event = dummy::GetEvent();
  EXPECT_THROW(ana->ProcessEvent(event), std::logic_error);
  auto reconstruction = new ReconstructionBase();
  reconstruction->Initialize(0);
  EXPECT_EQ(reconstruction->SelectEvent(event), true);
}

TEST(AnalysisTest, Clusterisation) {
  auto ana = std::make_unique<AnalysisBase>();
  auto track = dummy::GetTrack();
  auto clusters = ana->ClusterTrack(track);
  EXPECT_EQ(clusters.size(), 34);
  ana->setInvert(true);
  clusters = ana->ClusterTrack(track);
  EXPECT_EQ(clusters.size(), 3);
  ana->setClusterisation(nullptr);
  EXPECT_EXIT(ana->ClusterTrack(track), testing::ExitedWithCode(1), "ERROR");
}

//TEST(Graphics, drawEvent) {
//  auto ana = new AnalysisBase();
//  EXPECT_NO_THROW(ana->DrawSelection(GetEvent(), false));
//}