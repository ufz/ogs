append_source_files(SOURCES)
append_source_files(SOURCES BHE)
append_source_files(SOURCES BoundaryConditions)
append_source_files(SOURCES LocalAssemblers)

add_library(HeatTransportBHE ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatTransportBHE
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HeatTransportBHE PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
