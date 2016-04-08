include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Source files
APPEND_SOURCE_FILES(SOURCES)

add_subdirectory(GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib INTERFACE
    AssemblerLib
    MeshLib
    MeshGeoToolsLib
    NumLib # for shape matrices
    ${VTK_LIBRARIES}
    logog
)

ADD_VTK_DEPENDENCY(ProcessLib)

if(TARGET Eigen)
	add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(ProcessLib Boost)
endif()
