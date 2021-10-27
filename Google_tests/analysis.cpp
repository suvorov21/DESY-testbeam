//
// Created by SERGEY SUVOROV on 18/09/2021.
//
#include "gtest/gtest.h"
#include "AnalysisBase.hxx"

TEST(Initialisation, ParamFileRead) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->ReadParamFile(), true);
}

TEST(Initialisation, CLIread) {
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

TEST(Initialisation, outputWriteErr) {
  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->WriteOutput(), false);
}

TEST(InitialisationTest, baseProcessEvent) {
  auto ana = new AnalysisBase();
  auto event = new TEvent();
  EXPECT_THROW(ana->ProcessEvent(event), std::logic_error);
}

TEST(Graphics, drawEvent) {
  auto ana = new AnalysisBase();
  auto event = new TEvent();
  auto hit = new THit(10, 10, 200, 20);
  event->AddUsedHit(hit);
  EXPECT_NO_THROW(ana->DrawSelection(event));
}