append_source_files(SOURCES)

ogs_add_library(GroundwaterFlow ${SOURCES})

target_link_libraries(GroundwaterFlow PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
