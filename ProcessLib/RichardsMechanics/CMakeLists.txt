get_source_files(SOURCES)
append_source_files(SOURCES ConstitutiveRelations)

ogs_add_library(RichardsMechanics ${SOURCES})
target_link_libraries(RichardsMechanics PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(RichardsMechanics PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
