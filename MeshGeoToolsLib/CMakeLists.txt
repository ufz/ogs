# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)

# Create the library
add_library(MeshGeoToolsLib ${SOURCES_MeshGeoToolsLib})

target_link_libraries(MeshGeoToolsLib
    BaseLib
    MathLib
    MeshLib
    GeoLib
)

if(TARGET Boost)
    add_dependencies(MeshGeoToolsLib Boost)
endif()
