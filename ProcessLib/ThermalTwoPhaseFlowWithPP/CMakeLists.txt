APPEND_SOURCE_FILES(SOURCES)

add_library(ThermalTwoPhaseFlowWithPP ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowWithPP LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermalTwoPhaseFlowWithPP PUBLIC ProcessLib)

include(Tests.cmake)
