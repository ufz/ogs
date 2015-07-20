# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)

# Create the library
add_library(GeoLib STATIC ${SOURCES_GeoLib})

include_directories(
	.
)

target_link_libraries(GeoLib INTERFACE
	BaseLib
	MathLib
	logog
)

