APPEND_SOURCE_FILES(SOURCES)

add_library(GroundwaterFlow ${SOURCES})

target_link_libraries(GroundwaterFlow
    PUBLIC ProcessLib
)

include(Tests.cmake)
