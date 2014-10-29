
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/Legacy
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/VtkVis
)

IF(TARGET VtkVis)

	ADD_EXECUTABLE( generateStructuredQuadMesh generateStructuredQuadMesh.cpp )
	SET_TARGET_PROPERTIES( generateStructuredQuadMesh PROPERTIES FOLDER Utils)
	TARGET_LINK_LIBRARIES( generateStructuredQuadMesh
		MeshLib
		FileIO
		GeoLib
		MathLib
		VtkVis
		${CATALYST_LIBRARIES}
		${QT_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(generateStructuredQuadMesh)

	ADD_EXECUTABLE(  createMeshElemPropertiesFromASCRaster createMeshElemPropertiesFromASCRaster.cpp )
	SET_TARGET_PROPERTIES( createMeshElemPropertiesFromASCRaster PROPERTIES FOLDER Utils)
	TARGET_LINK_LIBRARIES( createMeshElemPropertiesFromASCRaster
		MeshLib
		FileIO
		VtkVis
		zlib
		${CATALYST_LIBRARIES}
		${QT_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(createMeshElemPropertiesFromASCRaster)

ENDIF() # VtkVis-target is existing


ADD_EXECUTABLE( generateStructuredMesh generateStructuredMesh.cpp )
TARGET_LINK_LIBRARIES( generateStructuredMesh
	BaseLib
	FileIO
	MathLib
	MeshLib
)
SET_TARGET_PROPERTIES(generateStructuredMesh PROPERTIES FOLDER Utilities)

