set(TOOLS computeSurfaceNodeIDsInPolygonalRegion constructMeshesFromGeometry
          identifySubdomains)
foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} MeshLib ApplicationsFileIO MeshGeoToolsLib)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)
