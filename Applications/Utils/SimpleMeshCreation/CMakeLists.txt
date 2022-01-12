set(TOOLS generateStructuredMesh createMeshElemPropertiesFromASCRaster)

foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} ApplicationsFileIO GitInfoLib MeshLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
