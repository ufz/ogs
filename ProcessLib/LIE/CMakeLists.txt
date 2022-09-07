add_subdirectory(Common)
get_source_files(SOURCES HydroMechanics)
append_source_files(SOURCES HydroMechanics/LocalAssembler)
append_source_files(SOURCES SmallDeformation)
append_source_files(SOURCES SmallDeformation/LocalAssembler)

ogs_add_library(LIE ${SOURCES})
target_link_libraries(LIE PUBLIC ProcessLib LIECommon PRIVATE ParameterLib)

target_precompile_headers(LIE PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(HydroMechanics/Tests.cmake)
    include(SmallDeformation/Tests.cmake)
endif()
