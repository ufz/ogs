# Source files
get_source_files(SOURCES_DataHolderLib)

# Library
ogs_add_library(DataHolderLib ${SOURCES_DataHolderLib})

target_link_libraries(
    DataHolderLib PUBLIC GeoLib MeshLib PRIVATE BaseLib MathLib spdlog::spdlog
)
