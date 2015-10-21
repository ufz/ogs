add_library(InSituLib
	VtkMappedMesh.h
	VtkMappedMesh.cpp
	VtkMappedMeshSource.h
	VtkMappedMeshSource.cpp
	VtkMappedPropertyVectorTemplate.h
	VtkMappedPropertyVectorTemplate-impl.h
	VtkMeshNodalCoordinatesTemplate.h
	VtkMeshNodalCoordinatesTemplate-impl.h
)

target_link_libraries(InSituLib INTERFACE MeshLib ${CATALYST_LIBRARIES})

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(InSituLib)
