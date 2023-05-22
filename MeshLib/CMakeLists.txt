get_source_files(SOURCES)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES NodePartitionedMesh.h)
endif()

append_source_files(SOURCES MeshSearch)
append_source_files(SOURCES Elements)
append_source_files(SOURCES IO)
append_source_files(SOURCES IO/Legacy)
append_source_files(SOURCES IO/VtkIO)
append_source_files(SOURCES Utils)
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES Utils/transformMeshToNodePartitionedMesh.cpp)
    list(REMOVE_ITEM SOURCES Utils/transformMeshToNodePartitionedMesh.h)
endif()
append_source_files(SOURCES IO/XDMF)
if(OGS_USE_PETSC)
    append_source_files(SOURCES IO/XDMF/mpi)
else()
    append_source_files(SOURCES IO/XDMF/posix)
endif()

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
           range-v3
           $<$<TARGET_EXISTS:VTK::IOParallelXML>:VTK::IOParallelXML>
           $<$<TARGET_EXISTS:VTK::ParallelMPI>:VTK::ParallelMPI>
           $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
           $<$<TARGET_EXISTS:petsc>:petsc>
    PRIVATE nlohmann_json::nlohmann_json
)

target_include_directories(MeshLib PUBLIC ${HDF5_INCLUDE_DIRS})

target_compile_definitions(
    MeshLib
    PUBLIC
        $<$<AND:$<BOOL:$ENV{VTK_USE_64BIT_IDS}>,$<NOT:$<BOOL:VTK_ADDED>>>:VTK_USE_64BIT_IDS>
)

target_precompile_headers(
    MeshLib
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    Elements/Element.h
    <Eigen/Core>
)
