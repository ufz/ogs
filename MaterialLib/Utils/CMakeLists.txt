get_source_files(SOURCES)

ogs_add_library(MaterialLib_Utils ${SOURCES})
target_link_libraries(MaterialLib_Utils PRIVATE BaseLib MeshLib)
