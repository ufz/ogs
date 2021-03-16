get_source_files(SOURCES)

ogs_add_library(ThermoRichardsMechanics ${SOURCES})
target_link_libraries(
    ThermoRichardsMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
