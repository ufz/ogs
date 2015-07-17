# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)

# Create the library
add_library(MeshGeoToolsLib STATIC ${SOURCES_MeshGeoToolsLib})

include_directories(
	.
	${CMAKE_CURRENT_SOURCE_DIR}/../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../GeoLib
)

target_link_libraries(MeshGeoToolsLib INTERFACE
	BaseLib
	MathLib
	MeshLib
	GeoLib
	logog
)

