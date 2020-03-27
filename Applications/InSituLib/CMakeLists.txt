# Source files
get_source_files(SOURCES)

# Library
ogs_add_library(InSituLib ${SOURCES})

target_link_libraries(PUBLIC BaseLib PRIVATE MeshLib)
