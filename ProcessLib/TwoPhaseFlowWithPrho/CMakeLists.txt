APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPrho ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TwoPhaseFlowWithPrho LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TwoPhaseFlowWithPrho PUBLIC ProcessLib)

include(Tests.cmake)
