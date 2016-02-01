add_executable(test_node_partitioned_mesh
	NodePartitionedMeshTester.cpp
)

target_link_libraries(test_node_partitioned_mesh
	MeshLib
	FileIO
	BaseLib
	logog
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
	EXECUTABLE_ARGS "mesh_3d ${Data_BINARY_DIR}/NodePartitionedMesh/ASCII"
	WRAPPER mpirun
	WRAPPER_ARGS "-np 3"
	TESTER diff
	DIFF_DATA mesh_3d_partition_0.msh mesh_3d_partition_1.msh mesh_3d_partition_2.msh
)

AddTest(
	NAME NodePartitionedMeshTestBinary
	PATH NodePartitionedMesh/Binary
	EXECUTABLE test_node_partitioned_mesh
	EXECUTABLE_ARGS "mesh_3d ${Data_BINARY_DIR}/NodePartitionedMesh/Binary"
	WRAPPER mpirun
	WRAPPER_ARGS "-np 3"
	TESTER diff
	DIFF_DATA mesh_3d_partition_0.msh mesh_3d_partition_1.msh mesh_3d_partition_2.msh
)

## FEM#PETSc solver tests
AddTest(
	NAME ParallelFEM_GroundWaterFlow2D
	PATH EllipticPETSc/
	EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj -gw_ksp_type bcgs -gw_pc_type bjacobi -gw_ksp_atol 1.e-10
	WRAPPER mpirun
	WRAPPER_ARGS "-np 3"
	TESTER diff
	DIFF_DATA
		quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_0.vtu
		quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_1.vtu
		quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_2.vtu
)

AddTest(
	NAME ParallelFEM_GroundWaterFlow3D_DirichletBC
	PATH EllipticPETSc/
	EXECUTABLE_ARGS cube_1e3.prj -gw_ksp_type bcgs -gw_pc_type bjacobi -gw_ksp_atol 1.e-10
	WRAPPER mpirun
	WRAPPER_ARGS "-np 3"
	TESTER diff
	DIFF_DATA
		cube_1e3_result_pcs_0_ts_1_0.vtu
		cube_1e3_result_pcs_0_ts_1_1.vtu
		cube_1e3_result_pcs_0_ts_1_2.vtu
)

AddTest(
	NAME ParallelFEM_GroundWaterFlow3D_NeumannBC
	PATH EllipticPETSc/
	EXECUTABLE_ARGS cube_1e3_neumann.prj -gw_ksp_type bcgs -gw_pc_type bjacobi -gw_ksp_atol 1.e-10
	WRAPPER mpirun
	WRAPPER_ARGS "-np 3"
	TESTER diff
	DIFF_DATA
		cube_1e3_neumann_result_pcs_0_ts_1_0.vtu
		cube_1e3_neumann_result_pcs_0_ts_1_1.vtu
		cube_1e3_neumann_result_pcs_0_ts_1_2.vtu
)
