include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_ASSEMBLERLIB)
set(SOURCES ${SOURCES_ASSEMBLERLIB})

# Create the library
add_library(AssemblerLib STATIC ${SOURCES})

target_link_libraries(AssemblerLib
	MeshLib
)

if(TARGET Eigen)
	add_dependencies(AssemblerLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(AssemblerLib Boost)
endif()
