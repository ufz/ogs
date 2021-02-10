get_source_files(SOURCES)

ogs_add_library(SmallDeformationNonlocal ${SOURCES})
target_link_libraries(SmallDeformationNonlocal
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
