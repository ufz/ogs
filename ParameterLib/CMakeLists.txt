get_source_files(SOURCES)
ogs_add_library(ParameterLib ${SOURCES})

target_link_libraries(ParameterLib
    PUBLIC MathLib
    PRIVATE BaseLib exprtk MeshLib spdlog::spdlog
)
