
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
