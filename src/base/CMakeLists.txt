message(STATUS "Base classes ...")

add_subdirectory(${PROJECT_SOURCE_DIR}/base/DataStorage)

file (GLOB_RECURSE BASE_SET_SOURCEFILES
    ${PROJECT_SOURCE_DIR}/base/AnalysisBase.cxx
    ${PROJECT_SOURCE_DIR}/base/ReconstructionBase.cxx
    ${PROJECT_SOURCE_DIR}/base/TrackFitter.cxx
    ${PROJECT_SOURCE_DIR}/base/TCluster.cxx
    ${PROJECT_SOURCE_DIR}/Reconstruction/*.cxx
    ${PROJECT_SOURCE_DIR}/utils/*.cxx
)

file (GLOB BASE_SET_HEADERFILES
    ${PROJECT_SOURCE_DIR}/base/AnalysisBase.hxx
    ${PROJECT_SOURCE_DIR}/base/ReconstructionBase.hxx
    ${PROJECT_SOURCE_DIR}/base/TrackFitter.hxx
    ${PROJECT_SOURCE_DIR}/base/TCluster.hxx
    ${PROJECT_SOURCE_DIR}/Reconstruction/*.hxx
    ${PROJECT_SOURCE_DIR}/utils/*.hxx
)

add_library(TBase ${BASE_SET_SOURCEFILES})
target_link_libraries(TBase TRawEvent)

add_custom_command(TARGET TBase POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                   ${CMAKE_CURRENT_BINARY_DIR}/DataStorage/libEventDict_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/../)
