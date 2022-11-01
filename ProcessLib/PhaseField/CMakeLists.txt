get_source_files(SOURCES)

ogs_add_library(PhaseField ${SOURCES})
target_link_libraries(PhaseField PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(PhaseField PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
