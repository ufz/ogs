# Source files
GET_SOURCE_FILES(SOURCES)

add_library(ProcessLib STATIC ${SOURCES})

target_link_libraries(ProcessLib INTERFACE
    AssemblerLib
    MeshLib
    MeshGeoToolsLib
    ${CATALYST_LIBRARIES}
    logog
)

ADD_VTK_DEPENDENCY(ProcessLib)

if(TARGET Eigen)
	add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(ProcessLib Boost)
endif()
