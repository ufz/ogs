# Source files
append_source_and_header_files(_sources _headers _private DIRECTORIES .)
append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    MeshEditing
    MeshGenerators
    MeshQuality
)

# Create the library
ogs_add_library(MeshToolsLib ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    MeshToolsLib
    PUBLIC
        GeoLib
        MeshLib
        range-v3
        VTK::IOXML
        VTK::IOLegacy
        $<$<AND:$<BOOL:${OGS_USE_PETSC}>,$<TARGET_EXISTS:VTK::IOParallelXML>>:VTK::IOParallelXML>
        $<$<AND:$<BOOL:${OGS_USE_PETSC}>,$<TARGET_EXISTS:VTK::ParallelMPI>>:VTK::ParallelMPI>
    PRIVATE BaseLib MathLib NumLib GitInfoLib
)

target_precompile_headers(
    MeshToolsLib
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    [["MeshLib/Mesh.h"]]
    [["MeshLib/Elements/Element.h"]]
)
