append_source_files(SOURCES)

add_library(HT ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS HT LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(HT PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
