INCLUDE_DIRECTORIES(
	.
	${CMAKE_SOURCE_DIR}/BaseLib/
	${CMAKE_SOURCE_DIR}/FileIO/
	${CMAKE_SOURCE_DIR}/MeshLib/
)

# Create the executable
ADD_EXECUTABLE( test_node_partitioned_mesh
	NodePartitionedMeshTester.cpp
)

TARGET_LINK_LIBRARIES ( test_node_partitioned_mesh
	MeshLib
	FileIO
	BaseLib
	logog
	${ADDITIONAL_LIBS}
	${BOOST_LIBRARIES}
)

IF (OGS_USE_PETSC)
	TARGET_LINK_LIBRARIES( test_node_partitioned_mesh  ${PETSC_LIBRARIES})
ENDIF (OGS_USE_PETSC)

IF (OGS_USE_MPI)
	TARGET_LINK_LIBRARIES( test_node_partitioned_mesh  ${MPI_CXX_LIBRARIES})
ENDIF ()

ExternalData_expand_arguments(
data
NodePartitionedMeshReader
DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/ASCII/Output/,REGEX:.*}
)

ExternalData_expand_arguments(
data
NodePartitionedMeshReader
DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/Binary/Output/,REGEX:.*}
)

set(FilePath "DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/ASCII/,REGEX:.*}")
set(MPITestParameters -np 3 "${PROJECT_BINARY_DIR}/bin/test_node_partitioned_mesh" "${FilePath}")
ExternalData_Add_Test(
    data
    NAME NodePartitionedMeshTestASCII
    COMMAND  "mpirun" ${MPITestParameters} 
)

set(FilePath "DATA{${ExternalData_SOURCE_ROOT}/NodePartitionedMesh/Binary/,REGEX:.*}")
set(MPITestParameters -np 3 "${PROJECT_BINARY_DIR}/bin/test_node_partitioned_mesh" "${FilePath}")
ExternalData_Add_Test(
    data
    NAME NodePartitionedMeshTestBinary
    COMMAND  "mpirun" ${MPITestParameters} 
)
