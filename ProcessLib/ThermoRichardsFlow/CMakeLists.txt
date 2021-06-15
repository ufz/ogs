get_source_files(SOURCES)

ogs_add_library(ThermoRichardsFlow ${SOURCES})
target_link_libraries(ThermoRichardsFlow PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
