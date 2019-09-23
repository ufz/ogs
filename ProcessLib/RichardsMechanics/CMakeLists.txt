append_source_files(SOURCES)

add_library(RichardsMechanics ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS RichardsMechanics
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(RichardsMechanics
                      PUBLIC ProcessLib
                      PRIVATE RichardsFlow ParameterLib)

include(Tests.cmake)
