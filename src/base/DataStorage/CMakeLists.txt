message(STATUS "TEvent dictionary ...")
#---Generate dictionaries
set (EVENT_DICT_HEADERFILES
        ${PROJECT_SOURCE_DIR}/src/base/DataStorage/*.hxx
)
file (GLOB EVENT_SOURCEFILES ${PROJECT_SOURCE_DIR}/src/base/DataStorage/*.cxx)

# Building a ROOT dictionary is an option
#set( EVENT_SOURCEFILES ${EVENT_SOURCEFILES} ${CMAKE_CURRENT_BINARY_DIR}/FullEventDict.cxx )
#ROOT_GENERATE_DICTIONARY(FullEventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${PROJECT_SOURCE_DIR}/src/base/DataStorage/linkDef.h)

#---Create a shared library with geneated dictionary
add_library(TEvent ${EVENT_SOURCEFILES})
target_link_libraries(TEvent ${ROOT_LIBRARIES} TRawEvent)