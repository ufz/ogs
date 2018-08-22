APPEND_SOURCE_FILES(SOURCES)

add_library(RichardsMechanics ${SOURCES})
target_link_libraries(RichardsMechanics PUBLIC ProcessLib RichardsFlow)

include(Tests.cmake)
