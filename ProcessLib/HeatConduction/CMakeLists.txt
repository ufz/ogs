get_source_files(SOURCES)

ogs_add_library(HeatConduction ${SOURCES})
target_link_libraries(HeatConduction PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
