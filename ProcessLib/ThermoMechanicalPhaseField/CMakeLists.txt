append_source_files(SOURCES)

add_library(ThermoMechanicalPhaseField ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermoMechanicalPhaseField
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermoMechanicalPhaseField
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

include(Tests.cmake)
