get_source_files(SOURCES)

ogs_add_library(SmallDeformation ${SOURCES})
target_link_libraries(SmallDeformation PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
