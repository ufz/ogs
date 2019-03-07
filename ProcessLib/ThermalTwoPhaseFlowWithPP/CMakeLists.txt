APPEND_SOURCE_FILES(SOURCES)

add_library(ThermalTwoPhaseFlowWithPP ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowWithPP LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                              RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(ThermalTwoPhaseFlowWithPP
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
