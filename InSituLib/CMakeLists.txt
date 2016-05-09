add_library(InSituLib STATIC
	VtkMappedMesh.h
	VtkMappedMesh.cpp
	VtkMappedMeshSource.h
	VtkMappedMeshSource.cpp
	VtkMappedPropertyVectorTemplate.h
	VtkMappedPropertyVectorTemplate-impl.h
	VtkMeshNodalCoordinatesTemplate.h
	VtkMeshNodalCoordinatesTemplate-impl.h
)

target_link_libraries(InSituLib MeshLib ${VTK_LIBRARIES})

ADD_VTK_DEPENDENCY(InSituLib)

if(TARGET Boost)
    add_dependencies(InSituLib Boost)
endif()
