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

  char* argv[2];
  argv[0] =  argv1;
  argv[1] =  argv2;

  auto ana = new AnalysisBase();
  EXPECT_EQ(ana->ReadCLI(argc, argv), true);
}