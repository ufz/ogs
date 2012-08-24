# Source files
GET_SOURCE_FILES(SOURCES_MESHLIB)
GET_SOURCE_FILES(SOURCES_ELEMENTS Elements)
GET_SOURCE_FILES(SOURCES_QUALITY MeshQuality)
SET ( SOURCES ${SOURCES_MESHLIB} ${SOURCES_ELEMENTS} ${SOURCES_QUALITY})

# Create the library
ADD_LIBRARY(MeshLib STATIC ${SOURCES})

include_directories(
	.
	../BaseLib
	../GeoLib
	../MathLib
	${CMAKE_SOURCE_DIR}/BaseLib/logog/include
)


target_link_libraries (MeshLib
	BaseLib
	GeoLib
	MathLib
	logog
)

