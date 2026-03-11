append_source_and_header_files(_sources _headers _private DIRECTORIES .)
ogs_add_library(ParameterLib ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    ParameterLib PUBLIC MathLib exprtk PRIVATE BaseLib MeshLib
)

target_precompile_headers(ParameterLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ParameterLib/Parameter.h"]] <Eigen/Core>)
