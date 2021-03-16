get_source_files(SOURCES)

ogs_add_library(ThermoMechanicalPhaseField ${SOURCES})
target_link_libraries(
    ThermoMechanicalPhaseField PUBLIC ProcessLib PRIVATE ParameterLib
)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
