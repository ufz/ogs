get_source_files(SOURCES)

ogs_add_library(RichardsComponentTransport ${SOURCES})
target_link_libraries(RichardsComponentTransport
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
