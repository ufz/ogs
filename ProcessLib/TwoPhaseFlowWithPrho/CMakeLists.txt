APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPrho ${SOURCES})
target_link_libraries(TwoPhaseFlowWithPrho PUBLIC ProcessLib)

include(Tests.cmake)
