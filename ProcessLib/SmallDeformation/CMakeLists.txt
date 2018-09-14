APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformation ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS SmallDeformation LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(SmallDeformation PUBLIC ProcessLib)

include(Tests.cmake)
