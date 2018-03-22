APPEND_SOURCE_FILES(SOURCES)

add_library(HydroMechanics ${SOURCES})
target_link_libraries(HydroMechanics PUBLIC ProcessLib)

include(Tests.cmake)
