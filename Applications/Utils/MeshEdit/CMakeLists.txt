set(TOOLS
    AddElementQuality
    AddFaultToVoxelGrid
    AddLayer
    appendLinesAlongPolyline
    checkMesh
    convertToLinearMesh
    CreateBoundaryConditionsAlongPolylines
    createLayeredMeshFromRasters
    createQuadraticMesh
    editMaterialID
    ExtractBoundary
    ExtractMaterials
    ExtractSurface
    Layers2Grid
    MapGeometryToMeshSurface
    MoveMesh
    moveMeshNodes
    NodeReordering
    queryMesh
    removeMeshElements
    ResetPropertiesInPolygonalRegion
    reviseMesh
    swapNodeCoordinateAxes
    Vtu2Grid)

if(OGS_BUILD_GUI)
    add_executable(RemoveGhostData
        RemoveGhostData.cpp
        ${PROJECT_SOURCE_DIR}/ThirdParty/paraview/vtkCleanUnstructuredGrid.cpp)
    target_include_directories(RemoveGhostData
        PRIVATE ${PROJECT_SOURCE_DIR}/ThirdParty/paraview)
    target_link_libraries(RemoveGhostData GitInfoLib MeshLib tclap)
    set_target_properties(RemoveGhostData PROPERTIES FOLDER Utilities)
    install(TARGETS RemoveGhostData RUNTIME DESTINATION bin COMPONENT Utilities)
endif()

foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} GitInfoLib MeshLib tclap)
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
