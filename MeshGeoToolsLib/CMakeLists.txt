# Source files
append_source_and_header_files(_sources _headers _private DIRECTORIES .)

# Create the library
ogs_add_library(MeshGeoToolsLib ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    MeshGeoToolsLib PUBLIC GeoLib MathLib
    PRIVATE BaseLib MeshToolsLib
            $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
)

target_precompile_headers(MeshGeoToolsLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] [["MeshLib/Mesh.h"]]
    [["MeshLib/Elements/Element.h"]])
