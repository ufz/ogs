get_source_files(SOURCES)

append_source_files(SOURCES ConstitutiveCommon)
append_source_files(SOURCES ConstitutiveOriginal)

ogs_add_library(ThermoRichardsMechanics ${SOURCES})
target_link_libraries(
    ThermoRichardsMechanics PUBLIC ProcessLib PRIVATE ParameterLib
)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
