get_source_files(SOURCES)

ogs_add_library(TES ${SOURCES})

target_link_libraries(TES PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(TES PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
