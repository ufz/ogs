#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_ASSEMBLERLIB)
SET ( SOURCES ${SOURCES_ASSEMBLERLIB})

# Create the library
ADD_LIBRARY( AssemblerLib STATIC ${SOURCES} )

INCLUDE_DIRECTORIES (
	.
	../BaseLib
	../GeoLib
	../MathLib
	../MeshLib
)

TARGET_LINK_LIBRARIES (AssemblerLib 
	${Boost_LIBRARIES}
	MeshLib
)

