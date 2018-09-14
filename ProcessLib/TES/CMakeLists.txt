APPEND_SOURCE_FILES(SOURCES)

add_library(TES ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TES LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TES PUBLIC ProcessLib)

include(Tests.cmake)
