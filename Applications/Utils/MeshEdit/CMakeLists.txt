
include_directories(
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/GmshIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

# Create executables
if(QT4_FOUND)

	add_executable(moveMeshNodes moveMeshNodes.cpp)
	target_link_libraries(moveMeshNodes FileIO)
	ADD_CATALYST_DEPENDENCY(moveMeshNodes)
	set_target_properties(moveMeshNodes PROPERTIES FOLDER Utilities)

	add_executable(MapGeometryToMeshSurface
		MapGeometryToMeshSurface.cpp )
	target_link_libraries(MapGeometryToMeshSurface FileIO MeshGeoToolsLib)
	ADD_CATALYST_DEPENDENCY(MapGeometryToMeshSurface)
	set_target_properties(MapGeometryToMeshSurface PROPERTIES FOLDER Utilities)

endif() # QT4_FOUND

add_executable(removeMeshElements removeMeshElements.cpp)
target_link_libraries(removeMeshElements FileIO)
ADD_CATALYST_DEPENDENCY(removeMeshElements)
set_target_properties(removeMeshElements PROPERTIES FOLDER Utilities)

add_executable(NodeReordering NodeReordering.cpp)
target_link_libraries(NodeReordering FileIO)
ADD_CATALYST_DEPENDENCY(NodeReordering)
set_target_properties(NodeReordering PROPERTIES FOLDER Utilities)

add_executable(MoveMesh MoveMesh.cpp)
target_link_libraries(MoveMesh FileIO)
ADD_CATALYST_DEPENDENCY(MoveMesh)
set_target_properties(MoveMesh PROPERTIES FOLDER Utilities)

add_executable(appendLinesAlongPolyline appendLinesAlongPolyline.cpp)
target_link_libraries(appendLinesAlongPolyline FileIO MeshGeoToolsLib)
ADD_CATALYST_DEPENDENCY(appendLinesAlongPolyline)
set_target_properties(appendLinesAlongPolyline PROPERTIES FOLDER Utilities)

add_executable(editMaterialID editMaterialID.cpp)
target_link_libraries(editMaterialID FileIO)
ADD_CATALYST_DEPENDENCY(editMaterialID)
set_target_properties(editMaterialID PROPERTIES FOLDER Utilities)

add_executable(checkMesh checkMesh.cpp)
target_link_libraries(checkMesh FileIO)
ADD_CATALYST_DEPENDENCY(checkMesh)
set_target_properties(checkMesh PROPERTIES FOLDER Utilities)

add_executable(reviseMesh reviseMesh.cpp)
target_link_libraries (reviseMesh FileIO)
ADD_CATALYST_DEPENDENCY(reviseMesh)
set_target_properties(reviseMesh PROPERTIES FOLDER Utilities)

add_executable(ResetPropertiesInPolygonalRegion
	ResetPropertiesInPolygonalRegion.cpp)
target_link_libraries(ResetPropertiesInPolygonalRegion FileIO)
set_target_properties(ResetPropertiesInPolygonalRegion
	PROPERTIES FOLDER Utilities)

add_executable(AddTopLayer AddTopLayer.cpp)
target_link_libraries(AddTopLayer FileIO)
set_target_properties(AddTopLayer PROPERTIES FOLDER Utilities)

add_executable(CreateBoundaryConditionsAlongPolylines
	CreateBoundaryConditionsAlongPolylines.cpp )
target_link_libraries(CreateBoundaryConditionsAlongPolylines FileIO MeshGeoToolsLib)
ADD_CATALYST_DEPENDENCY(CreateBoundaryConditionsAlongPolylines)
set_target_properties(CreateBoundaryConditionsAlongPolylines
	PROPERTIES FOLDER Utilities)

add_executable(queryMesh queryMesh.cpp)
target_link_libraries(queryMesh FileIO)
set_target_properties(queryMesh PROPERTIES FOLDER Utilities)

####################
### Installation ###
####################
install(TARGETS
	removeMeshElements
	NodeReordering
	MoveMesh
	appendLinesAlongPolyline
	editMaterialID
	checkMesh
	reviseMesh
	ResetPropertiesInPolygonalRegion
	AddTopLayer
	CreateBoundaryConditionsAlongPolylines
	queryMesh
	RUNTIME DESTINATION bin COMPONENT Utilities
)
if(QT4_FOUND)
	install(TARGETS
		moveMeshNodes
		MapGeometryToMeshSurface
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)
endif()
