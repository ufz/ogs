# Source files
GET_SOURCE_FILES(SOURCES_OGSLIB)
SET ( SOURCES ${SOURCES_OGSLIB})

# Create the library
ADD_LIBRARY(OgsLib STATIC ${SOURCES})

include_directories(
	.
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FemLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MeshLib
)


target_link_libraries (OgsLib
	BaseLib
	FemLib
	GeoLib
	MeshLib
)

