
include_directories(
	${CMAKE_SOURCE_DIR}/Utils/PartitionMesh
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MeshLib
)

add_executable(partmesh PartitionMesh.cpp MeshPartitioning.h MeshPartitioning.cpp)
set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh FileIO)
ADD_VTK_DEPENDENCY(partmesh)

####################
### Installation ###
####################
install(TARGETS partmesh RUNTIME DESTINATION bin COMPONENT ogs_partmesh)

cpack_add_component(ogs_partmesh
	DISPLAY_NAME "Mesh partitioning tool"
	DESCRIPTION "Mesh partitioning tool."
	GROUP Utilities
)
