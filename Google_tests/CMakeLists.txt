# 'Google_test' is the subproject name
project(Google_tests)

# 'lib' is the folder with Google Test sources
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_library(TFake FakeMaker.cpp)

# 'Google_Tests_run' is the target name
add_executable(Google_Tests_run analysis.cpp utils.cpp fitters.cpp)
target_link_libraries(Google_Tests_run TFake TSpatialResol)
target_link_libraries(Google_Tests_run gtest gtest_main)