APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BHE)
APPEND_SOURCE_FILES(SOURCES LocalAssemblers)

add_library(HeatTransportBHE ${SOURCES})
target_link_libraries(HeatTransportBHE PUBLIC ProcessLib)
if(BUILD_SHARED_LIBS)
    install(TARGETS HeatTransportBHE LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(Tests.cmake)
