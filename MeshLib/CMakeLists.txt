# Source files
GET_SOURCE_FILES(SOURCES_MESHLIB)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES_MESHLIB NodePartitionedMesh.h)
endif()

GET_SOURCE_FILES(SOURCES_ELEMENTS Elements)
GET_SOURCE_FILES(SOURCES_EDITING MeshEditing)
GET_SOURCE_FILES(SOURCES_GENERATORS MeshGenerators)
GET_SOURCE_FILES(SOURCES_SEARCH MeshSearch)
GET_SOURCE_FILES(SOURCES_IO IO)
GET_SOURCE_FILES(SOURCES_IO_LEGACY IO/Legacy)
GET_SOURCE_FILES(SOURCES_IO_VTKIO IO/VtkIO)
GET_SOURCE_FILES(SOURCES_QUALITY MeshQuality)
GET_SOURCE_FILES(SOURCES_VTK Vtk)
GET_SOURCE_FILES(SOURCES_UTILS Utils)

set(SOURCES ${SOURCES_MESHLIB} ${SOURCES_ELEMENTS} ${SOURCES_EDITING}
    ${SOURCES_GENERATORS} ${SOURCES_QUALITY} ${SOURCES_SEARCH}
    ${SOURCES_IO} ${SOURCES_IO_LEGACY}  ${SOURCES_IO_VTKIO} ${SOURCES_VTK}
    ${SOURCES_UTILS})

# It could be used for other MPI based DDC approach in future.
if(OGS_USE_PETSC)
    GET_SOURCE_FILES(SOURCES_MPI_IO IO/MPI_IO)
    set(SOURCES ${SOURCES} ${SOURCES_MPI_IO})
endif()

# Create the library
add_library(MeshLib ${SOURCES})

target_link_libraries(MeshLib
    BaseLib
    GeoLib
    MathLib
    ${VTK_LIBRARIES}
)

ADD_VTK_DEPENDENCY(MeshLib)

if(TARGET Eigen)
    add_dependencies(MeshLib Eigen)
endif()

if(TARGET Boost)
    add_dependencies(MeshLib Boost)
endif()
