get_source_files(SOURCES)

ogs_add_library(ThermalTwoPhaseFlowWithPP ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowWithPP
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
target_link_libraries(
    ThermalTwoPhaseFlowWithPP PUBLIC ProcessLib PRIVATE ParameterLib
)

target_precompile_headers(ThermalTwoPhaseFlowWithPP PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]]
    [["ProcessLib/Process.h"]] [["MaterialLib/MPL/Medium.h"]]
    [["MaterialLib/MPL/Property.h"]] <Eigen/Core>)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
