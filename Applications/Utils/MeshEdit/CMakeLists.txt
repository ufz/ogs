
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

	ADD_EXECUTABLE( moveMeshNodes moveMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( moveMeshNodes
		BaseLib
		FileIO
		MathLib
		MeshLib
		InSituLib
		${CATALYST_LIBRARIES}
		${QT_LIBRARIES}
	)

	SET_TARGET_PROPERTIES(moveMeshNodes
		PROPERTIES FOLDER Utilities)

ENDIF() # QT4_FOUND

ADD_EXECUTABLE( removeMeshElements removeMeshElements.cpp )
TARGET_LINK_LIBRARIES( removeMeshElements
	BaseLib
	FileIO
	MathLib
	MeshLib
	InSituLib
	${CATALYST_LIBRARIES}
)
SET_TARGET_PROPERTIES(removeMeshElements PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( DataExplorer5NodeReordering DataExplorer5NodeReordering.cpp )
TARGET_LINK_LIBRARIES( DataExplorer5NodeReordering
	FileIO
	MeshLib
	InSituLib
	${CATALYST_LIBRARIES}
)
SET_TARGET_PROPERTIES(DataExplorer5NodeReordering PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( MoveMesh MoveMesh.cpp )
TARGET_LINK_LIBRARIES( MoveMesh
	FileIO
	MeshLib
	InSituLib
	${CATALYST_LIBRARIES}
)
SET_TARGET_PROPERTIES(MoveMesh PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( appendLinesAlongPolyline appendLinesAlongPolyline.cpp )
TARGET_LINK_LIBRARIES( appendLinesAlongPolyline
	BaseLib
	FileIO
	MathLib
	MeshLib
	MeshGeoToolsLib
	InSituLib
	${CATALYST_LIBRARIES}
)

SET_TARGET_PROPERTIES(appendLinesAlongPolyline
	PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( editMaterialID editMaterialID.cpp )
TARGET_LINK_LIBRARIES( editMaterialID
	BaseLib
	FileIO
	MathLib
	MeshLib
	InSituLib
	${CATALYST_LIBRARIES}
)

SET_TARGET_PROPERTIES(editMaterialID
	PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( checkMesh checkMesh.cpp)
TARGET_LINK_LIBRARIES ( checkMesh
	BaseLib
	FileIO
	MathLib
	MeshLib
)

SET_TARGET_PROPERTIES(checkMesh
	PROPERTIES FOLDER Utilities)
