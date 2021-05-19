if(OGS_USE_MPI)
    return()
endif()

CPMAddPackage(
    NAME metis
    GITHUB_REPOSITORY scibuilder/metis
    GIT_TAG 982842a5ace9b3da2b2800817eb9e5fd3b42966b
    DOWNLOAD_ONLY YES
)
include(${PROJECT_SOURCE_DIR}/scripts/cmake/MetisSetup.cmake)

ogs_add_executable(
    partmesh PartitionMesh.cpp Metis.cpp NodeWiseMeshPartitioner.cpp
)
target_link_libraries(partmesh GitInfoLib MeshLib tclap)
add_dependencies(partmesh mpmetis)

ogs_add_executable(
    binaryToVTU BinaryToVTU.cpp NodeWiseMeshPartitioner.cpp
)
target_link_libraries(
    binaryToVTU
    GitInfoLib
    MeshLib
    tclap
    $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
)

install(TARGETS partmesh RUNTIME DESTINATION bin)
