APPEND_SOURCE_FILES(SOURCES)

add_library(HT ${SOURCES})
target_link_libraries(HT PUBLIC ProcessLib)

include(Tests.cmake)
