APPEND_SOURCE_FILES(SOURCES)

add_library(RichardsFlow ${SOURCES})
target_link_libraries(RichardsFlow PUBLIC ProcessLib)

include(Tests.cmake)
