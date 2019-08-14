set(TOOLS
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
    queryMesh
    removeMeshElements
    ResetPropertiesInPolygonalRegion
    reviseMesh
    swapNodeCoordinateAxes
    UnityPreprocessing)
foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} MeshLib)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)

target_link_libraries(appendLinesAlongPolyline MeshGeoToolsLib
                      ApplicationsFileIO)
target_link_libraries(CreateBoundaryConditionsAlongPolylines MeshGeoToolsLib
                      ApplicationsFileIO)
target_link_libraries(createLayeredMeshFromRasters ApplicationsFileIO)
target_link_libraries(MapGeometryToMeshSurface MeshGeoToolsLib)
target_link_libraries(ResetPropertiesInPolygonalRegion ApplicationsFileIO)
