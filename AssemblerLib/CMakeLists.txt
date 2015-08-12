#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_ASSEMBLERLIB)
set(SOURCES ${SOURCES_ASSEMBLERLIB})

# Create the library
add_library(AssemblerLib STATIC ${SOURCES})

target_link_libraries(AssemblerLib INTERFACE
	MeshLib
)

if(TARGET Boost)
	add_dependencies(AssemblerLib Boost)
endif()
