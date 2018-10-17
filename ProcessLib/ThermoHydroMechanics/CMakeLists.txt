APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoHydroMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermoHydroMechanics LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermoHydroMechanics PUBLIC ProcessLib)

include(Tests.cmake)
