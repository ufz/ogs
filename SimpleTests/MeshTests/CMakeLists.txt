
INCLUDE_DIRECTORIES(
        .
	${CMAKE_SOURCE_DIR}/Base/
	${CMAKE_SOURCE_DIR}/Base/logog/include
	${CMAKE_SOURCE_DIR}/FileIO/
	${CMAKE_SOURCE_DIR}/MathLib/
	${CMAKE_SOURCE_DIR}/MeshLib/
)

# Create the executable
ADD_EXECUTABLE( MeshRead
        MeshRead.cpp
        ${SOURCES}
        ${HEADERS}
)

TARGET_LINK_LIBRARIES ( MeshRead
	MeshLib
	FileIO	
	MathLib	
	Base
	GeoLib
)

