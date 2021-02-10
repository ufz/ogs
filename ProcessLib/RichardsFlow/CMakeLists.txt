get_source_files(SOURCES)

ogs_add_library(RichardsFlow ${SOURCES})
target_link_libraries(RichardsFlow PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
