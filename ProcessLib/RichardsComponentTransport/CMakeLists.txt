APPEND_SOURCE_FILES(SOURCES)

add_library(RichardsComponentTransport ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS RichardsComponentTransport LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(RichardsComponentTransport PUBLIC ProcessLib)

include(Tests.cmake)
