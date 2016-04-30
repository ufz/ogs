include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Source files
APPEND_SOURCE_FILES(SOURCES)

add_subdirectory(GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)

add_subdirectory(TES)
APPEND_SOURCE_FILES(SOURCES TES)

add_library(ProcessLib STATIC ${SOURCES})

target_link_libraries(ProcessLib INTERFACE
    AssemblerLib
    MaterialsLibAdsorption
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
