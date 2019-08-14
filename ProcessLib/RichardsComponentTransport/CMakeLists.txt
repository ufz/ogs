append_source_files(SOURCES)

add_library(RichardsComponentTransport ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS RichardsComponentTransport
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(RichardsComponentTransport
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

include(Tests.cmake)
