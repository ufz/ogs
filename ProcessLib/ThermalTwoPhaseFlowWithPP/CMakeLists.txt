ogs_add_process(ThermalTwoPhaseFlowWithPP)
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowWithPP
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
