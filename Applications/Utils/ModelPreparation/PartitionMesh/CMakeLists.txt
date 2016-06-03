###
add_executable(partmesh PartitionMesh.cpp MeshPartitioning.h MeshPartitioning.cpp)
set_target_properties(partmesh PROPERTIES FOLDER Utilities)
target_link_libraries(partmesh FileIO)

if(UTILS_PARTMESH_BUILD_WITH_METIS)
	include_directories(${METIS_PATH}/libmetis)
	include_directories(${METIS_PATH}/programs)

	set(METIS_SOURCES
	  ${METIS_PATH}/programs/mpmetis.c
	  ${METIS_PATH}/programs/cmdline_mpmetis.c
	  ${METIS_PATH}/programs/io.c
	  ${METIS_PATH}/programs/stat.c
	  metis_main.h
	)
	ADD_LIBRARY(METIS_Lib STATIC ${METIS_SOURCES} )

	target_link_libraries(partmesh METIS_Lib metis)
endif()

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
