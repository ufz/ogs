APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformation ${SOURCES})
target_link_libraries(SmallDeformation PUBLIC ProcessLib)

include(Tests.cmake)
