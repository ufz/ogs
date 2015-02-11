# Source files
GET_SOURCE_FILES(SOURCES_MESHLIB)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
	list(REMOVE_ITEM SOURCES_MESHLIB NodePartitionedMesh.h)
endif()

GET_SOURCE_FILES(SOURCES_ELEMENTS Elements)
GET_SOURCE_FILES(SOURCES_EDITING MeshEditing)
GET_SOURCE_FILES(SOURCES_GENERATORS MeshGenerators)
if(OGS_DONT_USE_VTK)
	list(REMOVE_ITEM SOURCES_GENERATORS MeshGenerators/VtkMeshConverter.h)
	list(REMOVE_ITEM SOURCES_GENERATORS MeshGenerators/VtkMeshConverter.cpp)
endif()

GET_SOURCE_FILES(SOURCES_QUALITY MeshQuality)
set(SOURCES ${SOURCES_MESHLIB} ${SOURCES_ELEMENTS} ${SOURCES_EDITING} ${SOURCES_GENERATORS} ${SOURCES_QUALITY})

# Create the library
add_library(MeshLib STATIC ${SOURCES})

target_link_libraries(MeshLib
	BaseLib
	GeoLib
	MathLib
	logog
)

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(MeshLib)
