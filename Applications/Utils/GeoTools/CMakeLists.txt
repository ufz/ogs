set(TOOLS MoveGeometry createRaster addDataToRaster generateGeometry RemoveUnusedPoints)
foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(
        ${tool} GeoLib GitInfoLib ApplicationsFileIO MeshToolsLib tclap
    )
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
