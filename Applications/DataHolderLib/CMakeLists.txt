# Source files
append_source_and_header_files(_sources _headers _private DIRECTORIES .)

# Library
ogs_add_library(DataHolderLib ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    DataHolderLib PUBLIC GeoLib MeshLib PRIVATE BaseLib MathLib
)
