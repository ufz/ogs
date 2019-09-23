# Source files
get_source_files(SOURCES)

# Create the library
add_library(MeshGeoToolsLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS MeshGeoToolsLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(MeshGeoToolsLib
                      PUBLIC GeoLib MathLib
                      PRIVATE BaseLib MeshLib logog)

if(OGS_USE_PCH)
    cotire(MeshGeoToolsLib)
endif()
