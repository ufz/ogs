APPEND_SOURCE_FILES(SOURCES)

add_library(LiquidFlow ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS LiquidFlow LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                               RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(LiquidFlow
    PUBLIC ProcessLib
    PRIVATE ParameterLib
)

include(Tests.cmake)
