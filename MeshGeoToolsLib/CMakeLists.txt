# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)

if(MSVC)
    include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}/../MathLib
        ${CMAKE_CURRENT_BINARY_DIR}/../MeshLib
        )
endif()

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
