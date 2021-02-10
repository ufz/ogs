get_source_files(SOURCES)

ogs_add_library(PhaseField ${SOURCES})
target_link_libraries(PhaseField PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
