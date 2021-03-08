get_source_files(SOURCES)

ogs_add_library(RichardsMechanics ${SOURCES})
target_link_libraries(RichardsMechanics PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
