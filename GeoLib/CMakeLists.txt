# Source files
GET_SOURCE_FILES(SOURCES_GeoLib)

# Create the library
ADD_LIBRARY(GeoLib STATIC ${SOURCES_GeoLib})


include_directories(
	.
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/MathLib
)


target_link_libraries (GeoLib
	BaseLib
	MathLib
)

