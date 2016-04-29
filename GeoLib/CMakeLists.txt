# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)
GET_SOURCE_FILES(SOURCES_IO IO)
GET_SOURCE_FILES(SOURCES_IO_Legacy IO/Legacy)

set(SOURCES ${SOURCES_GeoLib} ${SOURCES_IO} ${SOURCES_IO_Legacy})

# Create the library
add_library(GeoLib STATIC ${SOURCES}
	${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/tetgen/predicates.cxx
)

target_link_libraries(GeoLib
	BaseLib
	MathLib
)

if(TARGET Boost)
	add_dependencies(GeoLib Boost)
endif()
