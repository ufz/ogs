append_source_files(SOURCES)

add_library(ThermoMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermoMechanics LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermoMechanics PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
