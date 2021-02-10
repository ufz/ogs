get_source_files(SOURCES)

ogs_add_library(ThermalTwoPhaseFlowWithPP ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowWithPP
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermalTwoPhaseFlowWithPP
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
