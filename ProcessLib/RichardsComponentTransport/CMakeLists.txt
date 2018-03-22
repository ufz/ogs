APPEND_SOURCE_FILES(SOURCES)

add_library(RichardsComponentTransport ${SOURCES})
target_link_libraries(RichardsComponentTransport PUBLIC ProcessLib)

include(Tests.cmake)
