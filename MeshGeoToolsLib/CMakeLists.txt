# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)

# Create the library
add_library(MeshGeoToolsLib STATIC ${SOURCES_MeshGeoToolsLib})

target_link_libraries(MeshGeoToolsLib INTERFACE
	BaseLib
	MathLib
	MeshLib
	GeoLib
	FileIO
	logog
)

if(TARGET Boost)
	add_dependencies(MeshGeoToolsLib Boost)
endif()
