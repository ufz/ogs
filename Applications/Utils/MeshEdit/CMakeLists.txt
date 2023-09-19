set(TOOLS
    AddElementQuality
    AddFaultToVoxelGrid
    AddLayer
    appendLinesAlongPolyline
    checkMesh
    convertToLinearMesh
    CreateBoundaryConditionsAlongPolylines
    createLayeredMeshFromRasters
    createTetgenSmeshFromRasters
    createQuadraticMesh
    editMaterialID
    ExtractBoundary
    ExtractMaterials
    ExtractSurface
    ipDataToPointCloud
    Layers2Grid
    MapGeometryToMeshSurface
    MoveMesh
    MeshMapping
    NodeReordering
    queryMesh
    removeMeshElements
    ReorderMesh
    ResetPropertiesInPolygonalRegion
    reviseMesh
    swapNodeCoordinateAxes
    Vtu2Grid
)

add_subdirectory(PVTU2VTU)

ogs_add_executable(
    RemoveGhostData RemoveGhostData.cpp
    ${PROJECT_SOURCE_DIR}/ThirdParty/paraview/vtkCleanUnstructuredGrid.cpp
)
target_include_directories(
    RemoveGhostData PRIVATE ${PROJECT_SOURCE_DIR}/ThirdParty/paraview
)
target_link_libraries(
    RemoveGhostData GitInfoLib MeshLib MeshToolsLib tclap VTK::FiltersParallel
)
install(TARGETS RemoveGhostData RUNTIME DESTINATION bin)

foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} GitInfoLib MeshLib MeshToolsLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)

target_link_libraries(
    appendLinesAlongPolyline MeshGeoToolsLib ApplicationsFileIO
)
target_link_libraries(
    CreateBoundaryConditionsAlongPolylines MeshGeoToolsLib ApplicationsFileIO
)

target_link_libraries(createLayeredMeshFromRasters ApplicationsFileIO)
target_link_libraries(createTetgenSmeshFromRasters ApplicationsFileIO)
target_link_libraries(MapGeometryToMeshSurface MeshGeoToolsLib)
target_link_libraries(MeshMapping ApplicationsFileIO)
target_link_libraries(ResetPropertiesInPolygonalRegion ApplicationsFileIO)
target_link_libraries(ipDataToPointCloud NumLib)
