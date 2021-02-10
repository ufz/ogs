get_source_files(SOURCES)

ogs_add_library(HT ${SOURCES})
target_link_libraries(HT PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
