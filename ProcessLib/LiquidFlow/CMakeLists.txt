APPEND_SOURCE_FILES(SOURCES)

add_library(LiquidFlow ${SOURCES})
target_link_libraries(LiquidFlow PUBLIC ProcessLib)

include(Tests.cmake)
