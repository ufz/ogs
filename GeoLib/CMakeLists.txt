# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)

# Create the library
add_library(GeoLib ${SOURCES_GeoLib}
	${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/tetgen/predicates.cxx
)

target_link_libraries(GeoLib
	BaseLib
	MathLib
	logog
)

if(TARGET Boost)
	add_dependencies(GeoLib Boost)
endif()
