get_source_files(SOURCES)

ogs_add_library(ThermoHydroMechanics ${SOURCES})
target_link_libraries(
    ThermoHydroMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
