APPEND_SOURCE_FILES(SOURCES)

add_library(ThermalTwoPhaseFlowWithPP ${SOURCES})
target_link_libraries(ThermalTwoPhaseFlowWithPP PUBLIC ProcessLib)

include(Tests.cmake)
