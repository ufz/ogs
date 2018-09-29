add_executable(test_node_partitioned_mesh
    NodePartitionedMeshTester.cpp
)

target_link_libraries(test_node_partitioned_mesh
    MeshLib
    ${ADDITIONAL_LIBS}
    ${BOOST_LIBRARIES}
)

if(OGS_USE_PETSC)
    target_link_libraries(test_node_partitioned_mesh ${PETSC_LIBRARIES})
endif()

if(OGS_USE_MPI)
    target_link_libraries(test_node_partitioned_mesh ${MPI_CXX_LIBRARIES})
endif()

AddTest(
    NAME NodePartitionedMeshTestASCII
    PATH NodePartitionedMesh/ASCII
    EXECUTABLE test_node_partitioned_mesh
    EXECUTABLE_ARGS mesh_3d ${Data_BINARY_DIR}/NodePartitionedMesh/ASCII
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER diff
    DIFF_DATA mesh_3d_partition_0.msh mesh_3d_partition_1.msh mesh_3d_partition_2.msh
)

AddTest(
    NAME NodePartitionedMeshTestBinary
    PATH NodePartitionedMesh/Binary
    EXECUTABLE test_node_partitioned_mesh
    EXECUTABLE_ARGS mesh_3d ${Data_BINARY_DIR}/NodePartitionedMesh/Binary
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER diff
    DIFF_DATA mesh_3d_partition_0.msh mesh_3d_partition_1.msh mesh_3d_partition_2.msh
)
