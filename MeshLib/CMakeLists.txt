GET_SOURCE_FILES(SOURCES)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES NodePartitionedMesh.h)
endif()

APPEND_SOURCE_FILES(SOURCES MeshEditing)
APPEND_SOURCE_FILES(SOURCES MeshGenerators)
APPEND_SOURCE_FILES(SOURCES MeshSearch)
APPEND_SOURCE_FILES(SOURCES Elements)
APPEND_SOURCE_FILES(SOURCES IO)
APPEND_SOURCE_FILES(SOURCES IO/Legacy)
APPEND_SOURCE_FILES(SOURCES IO/VtkIO)
APPEND_SOURCE_FILES(SOURCES MeshQuality)
APPEND_SOURCE_FILES(SOURCES Vtk)

# It could be used for other MPI based DDC approach in future.
if(OGS_USE_PETSC)
    APPEND_SOURCE_FILES(SOURCES IO/MPI_IO)
endif()

# Create the library
add_library(MeshLib ${SOURCES})
if(OGS_USE_PCH)
    cotire(MeshLib)
endif()

target_link_libraries(MeshLib
    PUBLIC BaseLib GeoLib MathLib logog vtkIOXML
)

if(OGS_USE_MPI)
    target_link_libraries(MeshLib
        PUBLIC vtkIOParallelXML vtkParallelMPI
    )
endif()
