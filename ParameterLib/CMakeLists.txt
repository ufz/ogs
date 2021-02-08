get_source_files(SOURCES)
ogs_add_library(ParameterLib ${SOURCES})

include(GenerateExportHeader)
generate_export_header(ParameterLib)
target_include_directories(ParameterLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ParameterLib
    PUBLIC MathLib
    PRIVATE BaseLib exprtk MeshLib spdlog::spdlog
)
