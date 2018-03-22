APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoMechanics ${SOURCES})
target_link_libraries(ThermoMechanics PUBLIC ProcessLib)

include(Tests.cmake)
