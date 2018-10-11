if(COMPILER_IS_GCC AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "5.4")
    message(WARNING "gcc 5.4 segfaults when compiling partmesh!\n"
        "Disabled partmesh compilation!")
    return()
endif()
add_executable(partmesh PartitionMesh.cpp Metis.cpp NodeWiseMeshPartitioner.cpp)

set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh MeshLib)
add_dependencies(partmesh mpmetis)

####################
### Installation ###
####################
install(TARGETS partmesh RUNTIME DESTINATION bin COMPONENT ogs_partmesh)

cpack_add_component(ogs_partmesh
    DISPLAY_NAME "Mesh partitioning tool"
    DESCRIPTION "Mesh partitioning tool."
    GROUP Utilities
)
