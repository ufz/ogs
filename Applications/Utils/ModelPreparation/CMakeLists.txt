set(TOOLS ComputeNodeAreasFromSurfaceMesh convertVtkDataArrayToVtkDataArray
          createNeumannBc scaleProperty)
foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} MeshLib)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)

add_subdirectory(PartitionMesh)
