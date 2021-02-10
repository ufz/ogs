get_source_files(SOURCES)

ogs_add_library(TES ${SOURCES})

target_link_libraries(TES PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
