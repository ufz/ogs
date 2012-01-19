# Source files

GET_SOURCE_FILES(SOURCES)

# Create the library
ADD_LIBRARY( Base STATIC ${SOURCES})

SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        ../GeoLib
        ../MathLib
        .
)

