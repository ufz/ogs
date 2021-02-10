get_source_files(SOURCES)

ogs_add_library(HydroMechanics ${SOURCES})
target_link_libraries(HydroMechanics PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
