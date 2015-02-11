include_directories(
	.
	${CMAKE_SOURCE_DIR}/BaseLib/
	${CMAKE_SOURCE_DIR}/FileIO/
	${CMAKE_SOURCE_DIR}/MeshLib/
)

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

set(FilePath "DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/ASCII/,REGEX:.*}")
set(MPITestParameters -np 3 "${PROJECT_BINARY_DIR}/bin/test_node_partitioned_mesh" "${FilePath}/mesh_3d")
ExternalData_Add_Test(
	data
	NAME NodePartitionedMeshTestASCII
	COMMAND "mpirun" ${MPITestParameters}
)

set(FilePath "DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/Binary/,REGEX:.*}")
set(MPITestParameters -np 3 "${PROJECT_BINARY_DIR}/bin/test_node_partitioned_mesh" "${FilePath}/mesh_3d")
ExternalData_Add_Test(
	data
	NAME NodePartitionedMeshTestBinary
	COMMAND  "mpirun" ${MPITestParameters}
)
