APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPP ${SOURCES})
target_link_libraries(TwoPhaseFlowWithPP PUBLIC ProcessLib)

include(Tests.cmake)
