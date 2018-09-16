APPEND_SOURCE_FILES(SOURCES)

add_library(TwoPhaseFlowWithPP ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TwoPhaseFlowWithPP LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TwoPhaseFlowWithPP PUBLIC ProcessLib)

include(Tests.cmake)
