# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)
GET_SOURCE_FILES(SOURCES_SEARCH MeshSearch)

set(SOURCES ${SOURCES_MeshGeoToolsLib} ${SOURCES_SEARCH})

# Create the library
add_library(MeshGeoToolsLib ${SOURCES})

target_link_libraries(MeshGeoToolsLib
    BaseLib
    MathLib
    MeshLib
    GeoLib
)

if(TARGET Boost)
    add_dependencies(MeshGeoToolsLib Boost)
endif()
