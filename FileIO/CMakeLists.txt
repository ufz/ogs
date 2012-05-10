# Source files
GET_SOURCE_FILES(SOURCES_FILEIO)
SET ( SOURCES ${SOURCES_FILEIO})

# Create the library
ADD_LIBRARY(FileIO STATIC ${SOURCES})

include_directories(
	.
	../BaseLib
	../GeoLib
	../MathLib
	../MeshLib
)

target_link_libraries (FileIO
	GeoLib
	MeshLib
)

