message(STATUS "Base classes ...")

add_subdirectory(${PROJECT_SOURCE_DIR}/src/base/DataStorage)

file (GLOB_RECURSE BASE_SET_SOURCEFILES
    ${PROJECT_SOURCE_DIR}/src/base/*.cxx
    ${PROJECT_SOURCE_DIR}/src/Reconstruction/*.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/*.cxx
)

add_library(TBase ${BASE_SET_SOURCEFILES})
target_link_libraries(TBase TEvent)

#add_custom_command(TARGET TBase POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy
#                   ${CMAKE_CURRENT_BINARY_DIR}/DataStorage/libEventDict_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/../../)
