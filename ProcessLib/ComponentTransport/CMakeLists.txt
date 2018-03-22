APPEND_SOURCE_FILES(SOURCES)

add_library(ComponentTransport ${SOURCES})

target_link_libraries(ComponentTransport
    PUBLIC ProcessLib
)

include(Tests.cmake)
