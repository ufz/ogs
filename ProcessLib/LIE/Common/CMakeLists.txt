get_source_files(SOURCES)

ogs_add_library(LIECommon ${SOURCES})
target_link_libraries(LIECommon PUBLIC MeshLib)
