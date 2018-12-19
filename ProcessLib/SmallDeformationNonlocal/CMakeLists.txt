APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformationNonlocal ${SOURCES})
target_link_libraries(SmallDeformationNonlocal PUBLIC ProcessLib)

include(Tests.cmake)
