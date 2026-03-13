append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    MeshSearch
    Elements
    IO
    IO/Legacy
    IO/VtkIO
    IO/XDMF
    Utils
    Vtk
)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM _headers NodePartitionedMesh.h
         Utils/transformMeshToNodePartitionedMesh.h
    )
    list(REMOVE_ITEM _sources Utils/transformMeshToNodePartitionedMesh.cpp)
endif()
if(OGS_USE_PETSC)
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES IO/XDMF/mpi IO/MPI_IO
    )
else()
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES IO/XDMF/posix
    )
endif()

# Create the library
ogs_add_library(
    MeshLib ${_sources} PUBLIC_HEADERS ${_headers} PRIVATE_HEADERS ${_private}
)

target_link_libraries(
    MeshLib
    PUBLIC
        BaseLib
        GeoLib
        GitInfoLib
        ${HDF5_LIBRARIES}
        MathLib
        VTK::IOXML
        VTK::IOLegacy
        range-v3
        $<$<AND:$<BOOL:${OGS_USE_PETSC}>,$<TARGET_EXISTS:VTK::IOParallelXML>>:VTK::IOParallelXML>
        $<$<AND:$<BOOL:${OGS_USE_PETSC}>,$<TARGET_EXISTS:VTK::ParallelMPI>>:VTK::ParallelMPI>
        $<$<AND:$<BOOL:${OGS_USE_PETSC}>,$<TARGET_EXISTS:MPI::MPI_CXX>>:MPI::MPI_CXX>
        $<$<TARGET_EXISTS:petsc>:petsc>
    PRIVATE nlohmann_json::nlohmann_json
)

target_include_directories(MeshLib PUBLIC ${HDF5_INCLUDE_DIRS})

target_compile_definitions(
    MeshLib
    PUBLIC
        $<$<AND:$<BOOL:$ENV{VTK_USE_64BIT_IDS}>,$<NOT:$<BOOL:${VTK_ADDED}>>>:VTK_USE_64BIT_IDS>
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
