# Source files
GET_SOURCE_FILES(SOURCES)

# Create the library
add_library(MeshGeoToolsLib ${SOURCES})

target_link_libraries(MeshGeoToolsLib
    PUBLIC GeoLib MathLib
    PRIVATE BaseLib MeshLib logog
)

if(TARGET Boost)
    add_dependencies(MeshGeoToolsLib Boost)
endif()

if(OGS_USE_PCH)
    cotire(MeshGeoToolsLib)
endif()
