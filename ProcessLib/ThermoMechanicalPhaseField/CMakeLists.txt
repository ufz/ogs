APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoMechanicalPhaseField ${SOURCES})
target_link_libraries(ThermoMechanicalPhaseField PUBLIC ProcessLib)

include(Tests.cmake)
