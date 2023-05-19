# Check if cmake has the required version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)

# Set name of our project
project(DESY_testbeam VERSION 1.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(default_build_type "Debug")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_compile_options(-Wno-inconsistent-missing-override -Wno-nonportable-include-path)
add_link_options(-Wno-inconsistent-missing-override -Wno-nonportable-include-path)

if (NOT CMAKE_COVERAGE)
  message(STATUS "Switch off coverage flags")
  set(CMAKE_COVERAGE NO)
endif()
if (CMAKE_COVERAGE)
  message(STATUS "Switch on coverage flags")
  add_compile_options(-fprofile-arcs -ftest-coverage -Wno-inconsistent-missing-override)
  add_link_options(-fprofile-arcs -ftest-coverage -Wno-inconsistent-missing-override)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

if(ROOT_CXX_FLAGS MATCHES "-std=c\\+\\+17")
  set(CMAKE_CXX_STANDARD 17)
endif()

include_directories(external/cpp-generic-toolbox/include)
include_directories(external/simple-cpp-cmd-line-parser/include)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# Add hat_event package
set(PBUILDER OFF)
include_directories (
        ${ROOT_INCLUDE_DIR}
        # ${PROJECT_SOURCE_DIR}/external/hat_event
        ${PROJECT_SOURCE_DIR}/external/hat_event/src
)

add_subdirectory (external/hat_event)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/base/DataStorage) #---Event format
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/base/)            #---Base processors
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/Reconstruction/)  #---Reconstruction
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/utils/)           #---utils

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src//SpatialResol/)    #---SR analysis

add_subdirectory(${PROJECT_SOURCE_DIR}/src/base/)

#---Create  a main program using the library
add_library(TSpatialResol ${PROJECT_SOURCE_DIR}/src/SpatialResol/SpatialResolAna.cxx)
target_link_libraries(TSpatialResol TBase)
add_executable(SpatialResol.exe ${PROJECT_SOURCE_DIR}/app/SpatialResolutionApp.cxx)
target_link_libraries(SpatialResol.exe TSpatialResol)

# installation
install(TARGETS SpatialResol.exe TSpatialResol)

# Tests
add_subdirectory(Google_tests)
include_directories(Google_tests/lib/googletest/include)
