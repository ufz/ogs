append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(LIECommon ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(LIECommon PUBLIC MeshLib PRIVATE Boost::math)
