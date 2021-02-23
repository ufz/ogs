get_source_files(SOURCES)

ogs_add_library(ComponentTransport ${SOURCES})

target_link_libraries(ComponentTransport PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
