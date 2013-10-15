
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
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

ADD_EXECUTABLE( reverseMeshNodeOrdering reverseMeshNodeOrdering.cpp )
TARGET_LINK_LIBRARIES( reverseMeshNodeOrdering
	BaseLib
	FileIO
	MathLib
	MeshLib
)
SET_TARGET_PROPERTIES(reverseMeshNodeOrdering PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( removeMeshElements removeMeshElements.cpp )
TARGET_LINK_LIBRARIES( removeMeshElements
	BaseLib
	FileIO
	MathLib
	MeshLib
)
SET_TARGET_PROPERTIES(removeMeshElements PROPERTIES FOLDER Utilities)
