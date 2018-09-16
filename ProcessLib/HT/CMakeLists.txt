APPEND_SOURCE_FILES(SOURCES)

add_library(HT ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HT LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HT PUBLIC ProcessLib)

include(Tests.cmake)
