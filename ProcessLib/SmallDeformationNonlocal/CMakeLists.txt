get_source_files(SOURCES)

ogs_add_library(SmallDeformationNonlocal ${SOURCES})
target_link_libraries(
    SmallDeformationNonlocal PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(SmallDeformationNonlocal PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
