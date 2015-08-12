# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)

# Create the library
add_library(GeoLib STATIC ${SOURCES_GeoLib})

target_link_libraries(GeoLib INTERFACE
	BaseLib
	MathLib
	logog
)

if(TARGET Boost)
	add_dependencies(GeoLib Boost)
endif()
