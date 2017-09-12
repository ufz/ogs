add_executable(partmesh PartitionMesh.cpp NodeWiseMeshPartitioner.h NodeWiseMeshPartitioner.cpp)
set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh MeshLib)

####################
### Installation ###
####################
install(TARGETS partmesh RUNTIME DESTINATION bin COMPONENT ogs_partmesh)

cpack_add_component(ogs_partmesh
    DISPLAY_NAME "Mesh partitioning tool"
    DESCRIPTION "Mesh partitioning tool."
    GROUP Utilities
)
