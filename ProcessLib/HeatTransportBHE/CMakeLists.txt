APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BHE)
APPEND_SOURCE_FILES(SOURCES BoundaryConditions)
APPEND_SOURCE_FILES(SOURCES LocalAssemblers)

add_library(HeatTransportBHE ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatTransportBHE LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HeatTransportBHE
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
