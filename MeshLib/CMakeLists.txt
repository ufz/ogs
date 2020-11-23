get_source_files(SOURCES)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES NodePartitionedMesh.h)
endif()

append_source_files(SOURCES MeshEditing)
append_source_files(SOURCES MeshGenerators)
append_source_files(SOURCES MeshSearch)
append_source_files(SOURCES Elements)
append_source_files(SOURCES IO)
append_source_files(SOURCES IO/Legacy)
append_source_files(SOURCES IO/VtkIO)
append_source_files(SOURCES IO/XDMF)
append_source_files(SOURCES MeshQuality)
append_source_files(SOURCES Vtk)

# It could be used for other MPI based DDC approach in future.
if(OGS_USE_PETSC)
    append_source_files(SOURCES IO/MPI_IO)
endif()

# Create the library
ogs_add_library(MeshLib ${SOURCES})

target_link_libraries(
    MeshLib
    PUBLIC BaseLib GeoLib GitInfoLib MathLib ${VTK_LIBRARIES}
    PRIVATE spdlog::spdlog)

if(OGS_USE_MPI AND TARGET MPI::MPI_CXX)
    target_link_libraries(MeshLib PUBLIC MPI::MPI_CXX)
endif()
if(OGS_USE_PETSC)
    target_link_libraries(ogs PRIVATE ${PETSC_LIBRARIES})
endif()

if($ENV{VTK_USE_64BIT_IDS})
    target_compile_definitions(MeshLib PUBLIC VTK_USE_64BIT_IDS)
endif()
