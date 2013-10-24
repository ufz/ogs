# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)

# Create the library
ADD_LIBRARY(MeshGeoToolsLib STATIC ${SOURCES_MeshGeoToolsLib})

INCLUDE_DIRECTORIES(
	.
	${CMAKE_CURRENT_SOURCE_DIR}/../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../GeoLib
)

TARGET_LINK_LIBRARIES (MeshGeoToolsLib
	BaseLib
	MathLib
	MeshLib
	GeoLib
)

