get_source_files(SOURCES)

ogs_add_library(ThermoHydroMechanics ${SOURCES})
target_link_libraries(
    ThermoHydroMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(ThermoHydroMechanics PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
