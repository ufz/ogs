get_source_files(SOURCES)

ogs_add_library(TwoPhaseFlowWithPP ${SOURCES})
target_link_libraries(TwoPhaseFlowWithPP PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
