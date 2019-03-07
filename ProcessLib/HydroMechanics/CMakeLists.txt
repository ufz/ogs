APPEND_SOURCE_FILES(SOURCES)

add_library(HydroMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HydroMechanics LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(HydroMechanics
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
