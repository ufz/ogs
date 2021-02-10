get_source_files(SOURCES)

ogs_add_library(ThermoMechanics ${SOURCES})
target_link_libraries(ThermoMechanics PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
