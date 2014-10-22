INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_SOURCE_DIR}/../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib/Elements
)

ADD_LIBRARY(InSituLib
	VtkMappedMesh.h
	VtkMappedMesh.cpp
	VtkMappedMeshSource.h
	VtkMappedMeshSource.cpp
)

INCLUDE(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(InSituLib)
