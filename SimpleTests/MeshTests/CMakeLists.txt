IF (WIN32)
	SET ( ADDITIONAL_LIBS Winmm.lib )
ENDIF()

INCLUDE_DIRECTORIES(
	.
	${CMAKE_SOURCE_DIR}/BaseLib/
	${CMAKE_SOURCE_DIR}/FileIO/
	${CMAKE_SOURCE_DIR}/GeoLib/
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
	BaseLib
	GeoLib
	logog
	zlib
	${ADDITIONAL_LIBS}
	${BOOST_LIBRARIES}
)

# Create CollapseMeshNodes executable
ADD_EXECUTABLE( CollapseMeshNodes
	CollapseMeshNodes.cpp
	${SOURCES}
	${HEADERS}
)

TARGET_LINK_LIBRARIES ( CollapseMeshNodes
	MeshLib
	FileIO
	MathLib
	BaseLib
	GeoLib
	logog
	${BOOST_LIBRARIES}
)

# Create MeshSearchTest executable
ADD_EXECUTABLE( MeshSearchTest
	MeshSearchTest.cpp
	${SOURCES}
	${HEADERS}
)

TARGET_LINK_LIBRARIES ( MeshSearchTest
	MeshLib
	FileIO
	MathLib
	BaseLib
	GeoLib
	logog
	${ADDITIONAL_LIBS}
	${BOOST_LIBRARIES}
)
