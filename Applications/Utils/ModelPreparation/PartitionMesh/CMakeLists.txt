###
add_executable(partmesh PartitionMesh.cpp NodeWiseMeshPartitioner.h NodeWiseMeshPartitioner.cpp)
set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh MeshLib)

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

## Compile mpmetis
add_definitions(-DSVNINFO="${SVNREV}")
include_directories(${METIS_PATH}/libmetis)
include_directories(${METIS_PATH}/programs)
set(METIS_SOURCES
   ${METIS_PATH}/programs/mpmetis.c
   ${METIS_PATH}/programs/cmdline_mpmetis.c
   ${METIS_PATH}/programs/io.c
   ${METIS_PATH}/programs/stat.c
	)
add_executable(mpmetis ${METIS_SOURCES})

target_link_libraries(mpmetis metis)

install(TARGETS mpmetis RUNTIME DESTINATION bin COMPONENT ogs_metis)

cpack_add_component(ogs_metis
	DISPLAY_NAME "mpmetis"
	DESCRIPTION "METIS."
	GROUP Utilities
)
