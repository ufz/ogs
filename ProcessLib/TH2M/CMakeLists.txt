get_source_files(SOURCES)
append_source_files(SOURCES PhaseTransitionModels)
append_source_files(SOURCES ConstitutiveRelations)

ogs_add_library(TH2M ${SOURCES})
target_link_libraries(TH2M PUBLIC ProcessLib PRIVATE ParameterLib)

target_precompile_headers(TH2M PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
