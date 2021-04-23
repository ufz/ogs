get_source_files(SOURCES PhaseTransitionModels)

ogs_add_library(TH2M ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TH2M LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TH2M PUBLIC ProcessLib PRIVATE ParameterLib)
