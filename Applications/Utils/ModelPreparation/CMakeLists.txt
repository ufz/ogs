if(NOT OGS_USE_MPI)
    set(TOOLS ComputeNodeAreasFromSurfaceMesh convertVtkDataArrayToVtkDataArray
              createNeumannBc scaleProperty
    )
    foreach(tool ${TOOLS})
        ogs_add_executable(${tool} ${tool}.cpp)
        target_link_libraries(${tool} GitInfoLib MeshLib tclap)
    endforeach()
    install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
endif()

add_subdirectory(PartitionMesh)
