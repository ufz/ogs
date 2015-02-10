#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_ASSEMBLERLIB)
set ( SOURCES ${SOURCES_ASSEMBLERLIB})

# Create the library
add_library( AssemblerLib STATIC ${SOURCES} )

include_directories (
	.
	../BaseLib
	../GeoLib
	../MathLib
	../MeshLib
)

target_link_libraries (AssemblerLib 
	${Boost_LIBRARIES}
	MeshLib
)

