# Source files
get_source_files(SOURCES)

# Create the library
ogs_add_library(MeshGeoToolsLib ${SOURCES})

target_link_libraries(
    MeshGeoToolsLib PUBLIC GeoLib MathLib PRIVATE BaseLib MeshLib
)
