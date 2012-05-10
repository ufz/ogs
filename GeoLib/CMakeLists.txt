# Source files
GET_SOURCE_FILES(SOURCES_GEOLIB)

# Create the library
ADD_LIBRARY(GeoLib STATIC ${SOURCES_GEOLIB})


include_directories(
	.
	../BaseLib
	../MathLib
)


target_link_libraries (GeoLib
	BaseLib
	MathLib
)

