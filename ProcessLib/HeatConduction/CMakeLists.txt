append_source_files(SOURCES)

add_library(HeatConduction ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatConduction LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HeatConduction PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
