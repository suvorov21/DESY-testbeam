//
// Created by SERGEY SUVOROV on 18/09/2021.
//
#include "gtest/gtest.h"
#include "AnalysisBase.hxx"

auto GetDummyEvent() {
  auto event = new TEvent();
  auto hit = std::make_shared<THit>(10, 10, 200, 20);
  event->AddUsedHit(hit);
  return event;
}

TEST(InitialisationTest, ParamFileRead) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->ReadParamFile(), true);
}

TEST(InitialisationTest, ParamFileNotRead) {
  auto ana = new AnalysisBase();
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

  auto ana = new AnalysisBase();
  EXPECT_THROW(ana->ReadCLI(argc, argv_test1), std::logic_error);
  EXPECT_EQ(ana->ReadCLI(argc, argv_test2), true);
}

TEST(InitialisationTest, outputWriteErr) {
  auto ana = new AnalysisBase();
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
  auto ana = new AnalysisBase();
  auto event =GetDummyEvent();
  EXPECT_THROW(ana->ProcessEvent(event), std::logic_error);
  auto reconstruction = new ReconstructionBase();
  reconstruction->Initialize(0);
  EXPECT_EQ(reconstruction->SelectEvent(event), true);
}

//TEST(Graphics, drawEvent) {
//  auto ana = new AnalysisBase();
//  EXPECT_NO_THROW(ana->DrawSelection(GetDummyEvent(), false));
//}