APPEND_SOURCE_FILES(SOURCES)

add_library(ComponentTransport ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ComponentTransport LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

target_link_libraries(ComponentTransport
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
