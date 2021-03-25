# Check if cmake has the required version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)

# Set name of our project
project(DESY_testbeam VERSION 1.5)

set(default_build_type "Debug")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/base/DataStorage) #---Event format
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/base/)            #---Base processors
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/Reconstruction/)  #---Reconstruction
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/utils/)           #---utils

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/SpatialResol/)    #---SR analysis

add_subdirectory(${PROJECT_SOURCE_DIR}/base/)

#---Create  a main program using the library
add_executable(SpatialResol.exe ${PROJECT_SOURCE_DIR}/SpatialResol/SpatialResolAna.cxx ${PROJECT_SOURCE_DIR}/SpatialResol/SpatialResolAna.cxx)
target_link_libraries(SpatialResol.exe TRawEvent TBase)
