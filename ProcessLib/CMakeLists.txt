# Source files
GET_SOURCE_FILES(SOURCES)

add_library(ProcessLib STATIC ${SOURCES})

target_link_libraries(ProcessLib
    AssemblerLib
    MeshLib
    MeshGeoToolsLib
)
