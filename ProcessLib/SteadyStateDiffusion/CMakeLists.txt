append_source_files(SOURCES)

add_library(SteadyStateDiffusion ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS SteadyStateDiffusion LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(SteadyStateDiffusion PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
