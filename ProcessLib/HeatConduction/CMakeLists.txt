APPEND_SOURCE_FILES(SOURCES)

add_library(HeatConduction ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatConduction LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(HeatConduction
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
