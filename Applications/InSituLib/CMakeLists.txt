# Source files
GET_SOURCE_FILES(SOURCES)

# Library
add_library(InSituLib ${SOURCES})
if(OGS_USE_PCH)
    cotire(InSituLib)
endif()
