APPEND_SOURCE_FILES(SOURCES)

add_library(HeatConduction ${SOURCES})
target_link_libraries(HeatConduction PUBLIC ProcessLib)

include(Tests.cmake)
