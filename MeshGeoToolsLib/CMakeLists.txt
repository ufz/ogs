# Source files
get_source_files(SOURCES)

# Create the library
ogs_add_library(MeshGeoToolsLib ${SOURCES})

target_link_libraries(
    MeshGeoToolsLib PUBLIC GeoLib MathLib
    PRIVATE BaseLib MeshToolsLib
            $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
)

target_precompile_headers(MeshGeoToolsLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] [["MeshLib/Mesh.h"]]
    [["MeshLib/Elements/Element.h"]])
