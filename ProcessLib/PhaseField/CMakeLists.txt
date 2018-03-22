APPEND_SOURCE_FILES(SOURCES)

add_library(PhaseField ${SOURCES})
target_link_libraries(PhaseField PUBLIC ProcessLib)

include(Tests.cmake)
