CPMAddPackage(
    NAME metis
    GITHUB_REPOSITORY scibuilder/metis
    GIT_TAG 982842a5ace9b3da2b2800817eb9e5fd3b42966b
    DOWNLOAD_ONLY YES
)
include(${PROJECT_SOURCE_DIR}/scripts/cmake/MetisSetup.cmake)

add_executable(partmesh PartitionMesh.cpp Metis.cpp NodeWiseMeshPartitioner.cpp)
set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh GitInfoLib MeshLib tclap)
add_dependencies(partmesh mpmetis)
install(TARGETS partmesh RUNTIME DESTINATION bin)
