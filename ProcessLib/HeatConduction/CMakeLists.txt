APPEND_SOURCE_FILES(SOURCES)

add_library(HeatConduction ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatConduction LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HeatConduction PUBLIC ProcessLib)

include(Tests.cmake)
