# Source files
# GET_SOURCE_FILES(SOURCES_FILEIO)
SET( SOURCES
	Gmsh2GeoIO.cpp
	GMSInterface.cpp
	PetrelInterface.cpp
	readNonBlankLineFromInputStream.cpp
	Writer.cpp
)
GET_SOURCE_FILES(SOURCES_LEGACY Legacy)
GET_SOURCE_FILES(SOURCES_MESHIO MeshIO)
SET ( SOURCES ${SOURCES} ${SOURCES_LEGACY} ${SOURCES_MESHIO})

IF (QT4_FOUND)
	GET_SOURCE_FILES(SOURCES_XML XmlIO)
	SET ( SOURCES ${SOURCES} ${SOURCES_XML})
ENDIF (QT4_FOUND)

IF (Shapelib_FOUND)
	SET (SOURCES ${SOURCES} SHPInterface.h SHPInterface.cpp )
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
ENDIF (Shapelib_FOUND)

# Create the library
ADD_LIBRARY(FileIO STATIC ${SOURCES})

include_directories(
	.
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FemLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/OGS
	${CMAKE_BINARY_DIR}/BaseLib
)

target_link_libraries (FileIO
	GeoLib
	MeshLib
)

