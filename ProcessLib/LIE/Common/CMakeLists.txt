APPEND_SOURCE_FILES(SOURCES)

add_library(LIECommon ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS LIECommon LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                              RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
target_link_libraries(LIECommon PUBLIC MeshLib)
