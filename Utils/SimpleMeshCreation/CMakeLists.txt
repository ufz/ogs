
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/Legacy
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/VtkVis
)

INCLUDE( ${VTK_USE_FILE} )

# Create executables
IF(QT4_FOUND)
	ADD_EXECUTABLE( generateStructuredQuadMesh generateStructuredQuadMesh.cpp )
	SET_TARGET_PROPERTIES( generateStructuredQuadMesh PROPERTIES FOLDER Utils)
	TARGET_LINK_LIBRARIES( generateStructuredQuadMesh
		MeshLib
		FileIO
		VtkVis	
		vtkRendering
		${QT_LIBRARIES}
	)
ENDIF() # QT4_FOUND

