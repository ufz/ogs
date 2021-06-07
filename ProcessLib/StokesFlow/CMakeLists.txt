get_source_files(SOURCES)

ogs_add_library(StokesFlow ${SOURCES})

target_link_libraries(StokesFlow PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
