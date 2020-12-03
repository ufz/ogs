set(TOOLS
    AssignRasterDataToMesh
    computeSurfaceNodeIDsInPolygonalRegion
    constructMeshesFromGeometry
    createIntermediateRasters
    identifySubdomains
    IntegrateBoreholesIntoMesh
    Raster2Mesh
)

if(OGS_BUILD_GUI)
    # Utility requires that GUI is build because it uses Qt-based functionality
    # Also requires gmsh at runtime
    list(APPEND TOOLS VerticalSliceFromLayers)
endif()

foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL}
                          ApplicationsFileIO
                          GitInfoLib
                          MeshLib
                          MeshGeoToolsLib)
endforeach()

if(TARGET VerticalSliceFromLayers)
    target_link_libraries(VerticalSliceFromLayers Qt5::Xml Qt5::XmlPatterns)
endif()

set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)
