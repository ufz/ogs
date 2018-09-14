APPEND_SOURCE_FILES(SOURCES)

add_library(HydroMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HydroMechanics LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HydroMechanics PUBLIC ProcessLib)

include(Tests.cmake)
