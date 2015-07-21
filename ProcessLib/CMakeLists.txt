# Source files
GET_SOURCE_FILES(SOURCES)

add_library(ProcessLib STATIC ${SOURCES})

target_link_libraries(ProcessLib INTERFACE
    AssemblerLib
    MeshLib
    MeshGeoToolsLib
    logog
)

if(TARGET Eigen)
	add_dependencies(ProcessLib Eigen)
endif()
