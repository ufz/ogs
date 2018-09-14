APPEND_SOURCE_FILES(SOURCES)

add_library(RichardsFlow ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS RichardsFlow LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(RichardsFlow PUBLIC ProcessLib)

include(Tests.cmake)
