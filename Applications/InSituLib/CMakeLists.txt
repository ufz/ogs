# Source files
GET_SOURCE_FILES(SOURCES)

# Library
add_library(InSituLib ${SOURCES})
cotire(InSituLib)
