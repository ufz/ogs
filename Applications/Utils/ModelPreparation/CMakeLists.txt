set(TOOLS ComputeNodeAreasFromSurfaceMesh convertVtkDataArrayToVtkDataArray
          createNeumannBc scaleProperty
)
foreach(tool ${TOOLS})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} GitInfoLib MeshLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)

add_subdirectory(PartitionMesh)
