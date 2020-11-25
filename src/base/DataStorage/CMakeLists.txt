message(STATUS "TEvent dictionary ...")
#---Generate dictionaries
set (EVENT_DICT_HEADERFILES
    ${PROJECT_SOURCE_DIR}/base/DataStorage/TCluster.hxx
    ${PROJECT_SOURCE_DIR}/base/DataStorage/TEvent.hxx
    ${PROJECT_SOURCE_DIR}/base/DataStorage/TTrack.hxx
    ${PROJECT_SOURCE_DIR}/base/DataStorage/THit.hxx
)
file (GLOB EVENT_SOURCEFILES ${PROJECT_SOURCE_DIR}/base/DataStorage/*.cxx)
set( EVENT_SOURCEFILES ${EVENT_SOURCEFILES} ${CMAKE_CURRENT_BINARY_DIR}/EventDict.cxx )

ROOT_GENERATE_DICTIONARY(EventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${PROJECT_SOURCE_DIR}/base/DataStorage/linkDef.h)

#---Create a shared library with geneated dictionary
add_library(TEvent ${EVENT_SOURCEFILES})
target_link_libraries(TEvent ${ROOT_LIBRARIES})