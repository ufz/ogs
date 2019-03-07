# Source files
get_source_files(SOURCES)

add_library(ParameterLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ParameterLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

include(GenerateExportHeader)
generate_export_header(ParameterLib)
target_include_directories(ParameterLib
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(ParameterLib PRIVATE BaseLib MeshLib)

if(OGS_USE_PCH)
    cotire(ParameterLib)
endif()
