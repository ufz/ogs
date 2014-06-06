
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/GmshIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

# Create executables
IF(QT4_FOUND)

	ADD_EXECUTABLE( removeMeshNodes removeMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( removeMeshNodes
		BaseLib
		FileIO
		MeshLib
		${QT_LIBRARIES}
	)

	ADD_EXECUTABLE( moveMeshNodes moveMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( moveMeshNodes
		BaseLib
		FileIO
		MathLib
		MeshLib
		${QT_LIBRARIES}
	)

	SET_TARGET_PROPERTIES(moveMeshNodes removeMeshNodes
		PROPERTIES FOLDER Utilities)

ENDIF() # QT4_FOUND

ADD_EXECUTABLE( removeMeshElements removeMeshElements.cpp )
TARGET_LINK_LIBRARIES( removeMeshElements
	BaseLib
	FileIO
	MathLib
	MeshLib
)
SET_TARGET_PROPERTIES(removeMeshElements PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( DataExplorer5NodeReordering DataExplorer5NodeReordering.cpp )
TARGET_LINK_LIBRARIES( DataExplorer5NodeReordering
	FileIO
	MeshLib
)
SET_TARGET_PROPERTIES(DataExplorer5NodeReordering PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( MoveMesh MoveMesh.cpp )
TARGET_LINK_LIBRARIES( MoveMesh
	FileIO
	MeshLib
)
SET_TARGET_PROPERTIES(MoveMesh PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( appendLinesAlongPolyline appendLinesAlongPolyline.cpp )
TARGET_LINK_LIBRARIES( appendLinesAlongPolyline
	BaseLib
	FileIO
	MathLib
	MeshLib
	MeshGeoToolsLib
)

SET_TARGET_PROPERTIES(appendLinesAlongPolyline
	PROPERTIES FOLDER Utilities)
