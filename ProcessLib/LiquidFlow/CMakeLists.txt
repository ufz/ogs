APPEND_SOURCE_FILES(SOURCES)

add_library(LiquidFlow ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS LiquidFlow LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(LiquidFlow PUBLIC ProcessLib)

include(Tests.cmake)
