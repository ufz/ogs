APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPrho ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TwoPhaseFlowWithPrho LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(TwoPhaseFlowWithPrho
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
