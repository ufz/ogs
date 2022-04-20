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
append_source_files(SOURCES Utils)
append_source_files(SOURCES IO/XDMF)
if(OGS_USE_PETSC)
    append_source_files(SOURCES IO/XDMF/mpi)
else()
    append_source_files(SOURCES IO/XDMF/posix)
endif()

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
    PUBLIC BaseLib
           GeoLib
           GitInfoLib
           ${HDF5_LIBRARIES}
           MathLib
           VTK::IOXML
           VTK::IOLegacy
           $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
           $<$<TARGET_EXISTS:petsc>:petsc>
)

target_include_directories(MeshLib PUBLIC ${HDF5_INCLUDE_DIRS})

target_compile_definitions(
    MeshLib
    PUBLIC
        $<$<AND:$<BOOL:$ENV{VTK_USE_64BIT_IDS}>,$<NOT:$<BOOL:VTK_ADDED>>>:VTK_USE_64BIT_IDS>
)
