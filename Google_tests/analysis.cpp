//
// Created by SERGEY SUVOROV on 18/09/2021.
//
#include "AnalysisBase.hxx"
#include "FakeMaker.h"
#include "gtest/gtest.h"

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

    char *argv_test1[2];
    argv_test1[0] = argv1;
    argv_test1[1] = argv2;

    char *argv_test2[2];
    argv_test2[0] = argv1;
    argv_test2[1] = argv3;

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
    EXPECT_DEATH(ana->Initialize(), "No output file specified");
    ana->setOutputFile("output.root");
    EXPECT_DEATH(ana->Initialize(), "No input file specified");
}

TEST(InitialisationTest, baseProcessEvent) {
    auto ana = std::make_unique<AnalysisBase>();
    auto event = fake::GetEvent();
    EXPECT_THROW(ana->ProcessEvent(event), std::logic_error);
    auto reconstruction = new ReconstructionBase();
    reconstruction->Initialize(0);
    EXPECT_EQ(reconstruction->ReconstructEvent(event), true);
}

TEST(AnalysisTest, Clusterisation) {
    Clustering cluster;
    auto track = fake::GetTrack();
    auto clusters = cluster.ClusterTrack(track);
    EXPECT_EQ(clusters.size(), 36);
    cluster = Clustering(clusterType::kRowColumn, true);
    clusters = cluster.ClusterTrack(track);
    EXPECT_EQ(clusters.size(), 3);
}

//TEST(Graphics, drawEvent) {
//  auto ana = new AnalysisBase();
//  EXPECT_NO_THROW(ana->DrawSelection(GetEvent(), false));
//}