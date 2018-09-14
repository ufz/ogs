APPEND_SOURCE_FILES(SOURCES)

add_library(PhaseField ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS PhaseField LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(PhaseField PUBLIC ProcessLib)

include(Tests.cmake)
