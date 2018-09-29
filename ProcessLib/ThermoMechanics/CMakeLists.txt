APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermoMechanics LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermoMechanics PUBLIC ProcessLib)

include(Tests.cmake)
