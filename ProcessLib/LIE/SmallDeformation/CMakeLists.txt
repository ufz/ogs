get_source_files(SOURCES)
append_source_files(SOURCES LocalAssembler)

ogs_add_library(LIE_M ${SOURCES})
target_link_libraries(LIE_M PUBLIC ProcessLib LIECommon PRIVATE ParameterLib)

target_precompile_headers(LIE_M PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
