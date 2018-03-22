APPEND_SOURCE_FILES(SOURCES)

add_library(TES ${SOURCES})
target_link_libraries(TES PUBLIC ProcessLib)

include(Tests.cmake)
