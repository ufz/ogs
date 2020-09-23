set(TOOLS
    AssignRasterDataToMesh
    computeSurfaceNodeIDsInPolygonalRegion
    constructMeshesFromGeometry
    createIntermediateRasters
    identifySubdomains
    IntegrateBoreholesIntoMesh
    Raster2Mesh
    VerticalSliceFromLayers
)
foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL}
                          ApplicationsFileIO
                          GitInfoLib
                          MeshLib
                          MeshGeoToolsLib)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)
