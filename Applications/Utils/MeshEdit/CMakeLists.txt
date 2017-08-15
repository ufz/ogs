
add_executable(moveMeshNodes moveMeshNodes.cpp)
target_link_libraries(moveMeshNodes MeshLib)
ADD_VTK_DEPENDENCY(moveMeshNodes)
set_target_properties(moveMeshNodes PROPERTIES FOLDER Utilities)

add_executable(MapGeometryToMeshSurface
    MapGeometryToMeshSurface.cpp )
target_link_libraries(MapGeometryToMeshSurface MeshGeoToolsLib MeshLib)
ADD_VTK_DEPENDENCY(MapGeometryToMeshSurface)
set_target_properties(MapGeometryToMeshSurface PROPERTIES FOLDER Utilities)

add_executable(removeMeshElements removeMeshElements.cpp)
target_link_libraries(removeMeshElements MeshLib)
ADD_VTK_DEPENDENCY(removeMeshElements)
set_target_properties(removeMeshElements PROPERTIES FOLDER Utilities)

add_executable(NodeReordering NodeReordering.cpp)
target_link_libraries(NodeReordering MeshLib)
ADD_VTK_DEPENDENCY(NodeReordering)
set_target_properties(NodeReordering PROPERTIES FOLDER Utilities)

add_executable(MoveMesh MoveMesh.cpp)
target_link_libraries(MoveMesh MeshLib)
ADD_VTK_DEPENDENCY(MoveMesh)
set_target_properties(MoveMesh PROPERTIES FOLDER Utilities)

add_executable(appendLinesAlongPolyline appendLinesAlongPolyline.cpp)
target_link_libraries(appendLinesAlongPolyline
    MeshGeoToolsLib
    MeshLib
    ApplicationsFileIO
)
ADD_VTK_DEPENDENCY(appendLinesAlongPolyline)
set_target_properties(appendLinesAlongPolyline PROPERTIES FOLDER Utilities)

add_executable(editMaterialID editMaterialID.cpp)
target_link_libraries(editMaterialID MeshLib)
ADD_VTK_DEPENDENCY(editMaterialID)
set_target_properties(editMaterialID PROPERTIES FOLDER Utilities)

add_executable(checkMesh checkMesh.cpp)
target_link_libraries(checkMesh MeshLib)
ADD_VTK_DEPENDENCY(checkMesh)
set_target_properties(checkMesh PROPERTIES FOLDER Utilities)

add_executable(reviseMesh reviseMesh.cpp)
target_link_libraries (reviseMesh MeshLib)
ADD_VTK_DEPENDENCY(reviseMesh)
set_target_properties(reviseMesh PROPERTIES FOLDER Utilities)

add_executable(ResetPropertiesInPolygonalRegion
    ResetPropertiesInPolygonalRegion.cpp)
target_link_libraries(ResetPropertiesInPolygonalRegion
    MeshLib
    ApplicationsFileIO
)
set_target_properties(ResetPropertiesInPolygonalRegion
    PROPERTIES FOLDER Utilities)

add_executable(AddTopLayer AddTopLayer.cpp)
target_link_libraries(AddTopLayer MeshLib)
set_target_properties(AddTopLayer PROPERTIES FOLDER Utilities)

add_executable(createLayeredMeshFromRasters createLayeredMeshFromRasters.cpp)
target_link_libraries(createLayeredMeshFromRasters
    MeshLib
    ApplicationsFileIO
)
set_target_properties(createLayeredMeshFromRasters PROPERTIES FOLDER Utilities)

add_executable(CreateBoundaryConditionsAlongPolylines
    CreateBoundaryConditionsAlongPolylines.cpp )
target_link_libraries(CreateBoundaryConditionsAlongPolylines
    MeshGeoToolsLib
    MeshLib
    ApplicationsFileIO
)
ADD_VTK_DEPENDENCY(CreateBoundaryConditionsAlongPolylines)
set_target_properties(CreateBoundaryConditionsAlongPolylines
    PROPERTIES FOLDER Utilities)

add_executable(queryMesh queryMesh.cpp)
target_link_libraries(queryMesh MeshLib)
set_target_properties(queryMesh PROPERTIES FOLDER Utilities)

add_executable(ExtractSurface ExtractSurface.cpp)
target_link_libraries(ExtractSurface MeshLib)
set_target_properties(ExtractSurface PROPERTIES FOLDER Utilities)

add_executable(swapNodeCoordinateAxes swapNodeCoordinateAxes.cpp)
target_link_libraries(swapNodeCoordinateAxes MeshLib)
set_target_properties(swapNodeCoordinateAxes PROPERTIES FOLDER Utilities)

add_executable(createQuadraticMesh createQuadraticMesh.cpp )
target_link_libraries(createQuadraticMesh MeshLib)
set_target_properties(createQuadraticMesh PROPERTIES FOLDER Utilities)

add_executable(convertToLinearMesh convertToLinearMesh.cpp )
target_link_libraries(convertToLinearMesh MeshLib)
set_target_properties(convertToLinearMesh PROPERTIES FOLDER Utilities)

add_executable(UnityPreprocessing UnityPreprocessing.cpp)
target_link_libraries(UnityPreprocessing MeshLib)
set_target_properties(UnityPreprocessing PROPERTIES FOLDER Utilities)

####################
### Installation ###
####################
install(TARGETS
    AddTopLayer
    appendLinesAlongPolyline
    checkMesh
    convertToLinearMesh
    CreateBoundaryConditionsAlongPolylines
    createLayeredMeshFromRasters
    createQuadraticMesh
    editMaterialID
    ExtractSurface
    MapGeometryToMeshSurface
    MoveMesh
    moveMeshNodes
    NodeReordering
    removeMeshElements
    ResetPropertiesInPolygonalRegion
    reviseMesh
    queryMesh
    swapNodeCoordinateAxes
    UnityPreprocessing
    RUNTIME DESTINATION bin COMPONENT Utilities
)
