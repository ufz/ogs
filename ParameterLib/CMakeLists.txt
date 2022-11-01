get_source_files(SOURCES)
ogs_add_library(ParameterLib ${SOURCES})

target_link_libraries(
    ParameterLib PUBLIC MathLib PRIVATE BaseLib exprtk MeshLib
)

target_precompile_headers(ParameterLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ParameterLib/Parameter.h"]] <Eigen/Core>)
