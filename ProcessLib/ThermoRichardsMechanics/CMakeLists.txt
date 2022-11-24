get_source_files(SOURCES)

append_source_files(SOURCES ConstitutiveCommon)
append_source_files(SOURCES ConstitutiveStress_StrainTemperature)
if (OGS_USE_MFRONT)
    append_source_files(SOURCES
                        ConstitutiveStressSaturation_StrainPressureTemperature)
endif()

ogs_add_library(ThermoRichardsMechanics ${SOURCES})
target_link_libraries(
    ThermoRichardsMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(ThermoRichardsMechanics PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
