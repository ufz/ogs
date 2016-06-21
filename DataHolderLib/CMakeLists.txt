# Source files
GET_SOURCE_FILES(SOURCES_DataHolderLib)

if(MSVC)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/../MathLib)
endif()

# Library
add_library(DataHolderLib ${SOURCES_DataHolderLib})

target_link_libraries(DataHolderLib
    GeoLib
    MeshLib
    logog
)