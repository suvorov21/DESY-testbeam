message(STATUS "TEvent dictionary ...")
#---Generate dictionaries
set (EVENT_DICT_HEADERFILES
    ${PROJECT_SOURCE_DIR}/src/base/DataStorage/TRawEvent.hxx
    ${PROJECT_SOURCE_DIR}/src/base/DataStorage/THit.hxx
)
file (GLOB EVENT_SOURCEFILES ${PROJECT_SOURCE_DIR}/src/base/DataStorage/*.cxx)
set( EVENT_SOURCEFILES ${EVENT_SOURCEFILES} ${CMAKE_CURRENT_BINARY_DIR}/EventDict.cxx )

ROOT_GENERATE_DICTIONARY(EventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${PROJECT_SOURCE_DIR}/src/base/DataStorage/linkDef.h)

#---Create a shared library with geneated dictionary
add_library(TRawEvent ${EVENT_SOURCEFILES})
target_link_libraries(TRawEvent ${ROOT_LIBRARIES})