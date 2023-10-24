if(OGS_USE_MPI)
    ogs_add_executable(
        binaryToPVTU BinaryToPVTU.cpp NodeWiseMeshPartitioner.cpp
    )
    target_link_libraries(
           binaryToPVTU GitInfoLib ProcessLib MeshLib MeshToolsLib tclap
        $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX> VTK::ParallelMPI
    )
    install(TARGETS binaryToPVTU RUNTIME DESTINATION bin)
endif()

ogs_add_executable(
    partmesh PartitionMesh.cpp Metis.cpp NodeWiseMeshPartitioner.cpp
)
target_link_libraries(partmesh GitInfoLib MeshLib MeshToolsLib tclap NumLib)
if(TARGET mpmetis)
    add_dependencies(partmesh mpmetis)
endif()
install(TARGETS partmesh RUNTIME DESTINATION bin)
