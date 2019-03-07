APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformation ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS SmallDeformation LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(SmallDeformation
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
