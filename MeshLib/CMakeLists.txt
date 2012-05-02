# Source files
GET_SOURCE_FILES(SOURCES_MESHLIB)
SET ( SOURCES ${SOURCES_MESHLIB})

GET_SOURCE_FILES(SOURCES_ELEMENTS Elements)
SET ( SOURCES ${SOURCES} ${SOURCES_ELEMENTS})

# Create the library
ADD_LIBRARY(MeshLib STATIC ${SOURCES})

include_directories(
	.
	../Base
	../GeoLib
	../MathLib
)


target_link_libraries (MeshLib
	Base
	GeoLib
	MathLib
)

