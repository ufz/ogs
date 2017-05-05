# Source files
GET_SOURCE_FILES(SOURCES_DataHolderLib)

# Library
add_library(DataHolderLib ${SOURCES_DataHolderLib})

target_link_libraries(DataHolderLib
    PUBLIC GeoLib MeshLib
    PRIVATE BaseLib logog
)

if(OGS_USE_PCH)
    cotire(DataHolderLib)
endif()
