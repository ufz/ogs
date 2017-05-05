# Source files
GET_SOURCE_FILES(SOURCES)

# Library
add_library(InSituLib ${SOURCES})

target_link_libraries(PUBLIC BaseLib PRIVATE MeshLib)

if(OGS_USE_PCH)
    cotire(InSituLib)
endif()
