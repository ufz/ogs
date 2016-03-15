include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl-dynamic.h.in"
		"${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl.h" COPYONLY
	)
else()
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl-fixed.h.in"
		"${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl.h" COPYONLY
	)
endif()

# Source files
APPEND_SOURCE_FILES(SOURCES)

add_subdirectory(GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)

add_library(ProcessLib STATIC ${SOURCES})

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
