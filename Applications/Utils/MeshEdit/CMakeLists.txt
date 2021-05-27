if(OGS_USE_MPI)
    return()
endif()

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
    Vtu2Grid
)

if(OGS_BUILD_GUI)
    ogs_add_executable(
        RemoveGhostData RemoveGhostData.cpp
        ${PROJECT_SOURCE_DIR}/ThirdParty/paraview/vtkCleanUnstructuredGrid.cpp
    )
    target_include_directories(
        RemoveGhostData PRIVATE ${PROJECT_SOURCE_DIR}/ThirdParty/paraview
    )
    target_link_libraries(RemoveGhostData GitInfoLib MeshLib tclap)
    install(TARGETS RemoveGhostData RUNTIME DESTINATION bin)
endif()

foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} GitInfoLib MeshLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)

target_link_libraries(
    appendLinesAlongPolyline MeshGeoToolsLib ApplicationsFileIO
)
target_link_libraries(
    CreateBoundaryConditionsAlongPolylines MeshGeoToolsLib ApplicationsFileIO
)
target_link_libraries(createLayeredMeshFromRasters ApplicationsFileIO)
target_link_libraries(MapGeometryToMeshSurface MeshGeoToolsLib)
target_link_libraries(ResetPropertiesInPolygonalRegion ApplicationsFileIO)
