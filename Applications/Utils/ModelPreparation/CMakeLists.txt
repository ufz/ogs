add_subdirectory(PartitionMesh)
add_subdirectory(CreateAnchors)

set(TOOLS ComputeNodeAreasFromSurfaceMesh convertVtkDataArrayToVtkDataArray
          createNeumannBc scaleProperty
)
foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} GitInfoLib MeshLib MeshToolsLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
