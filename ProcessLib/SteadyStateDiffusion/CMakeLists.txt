append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(SteadyStateDiffusion ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    SteadyStateDiffusion PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(SteadyStateDiffusion PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
