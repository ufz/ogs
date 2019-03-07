APPEND_SOURCE_FILES(SOURCES)

add_library(TES ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TES LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(TES
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
