append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(MaterialLib_Utils ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(MaterialLib_Utils PUBLIC BaseLib PRIVATE MeshLib)
