append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(
    MaterialLib_FractureModels ${_sources} GENERATE_EXPORT_HEADER
    PUBLIC_HEADERS ${_headers}
)
target_link_libraries(
    MaterialLib_FractureModels PUBLIC MathLib PRIVATE ParameterLib NumLib
                                                      Boost::math
)
