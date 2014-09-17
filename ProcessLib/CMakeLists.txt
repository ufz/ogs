# Source files
GET_SOURCE_FILES(SOURCES)

ADD_LIBRARY(ProcessLib STATIC ${SOURCES})

TARGET_LINK_LIBRARIES( ProcessLib
    AssemblerLib
    MeshLib
)
