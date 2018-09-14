APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoMechanicalPhaseField ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermoMechanicalPhaseField LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermoMechanicalPhaseField PUBLIC ProcessLib)

include(Tests.cmake)
