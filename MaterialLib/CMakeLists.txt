# Source files
GET_SOURCE_FILES(SOURCES)
GET_SOURCE_FILES(SOURCES_ADSORPTION Adsorption)
set(SOURCES ${SOURCES} ${SOURCES_ADSORPTION})

APPEND_SOURCE_FILES(SOURCES SolidModels)

add_library(MaterialLib ${SOURCES} )
target_link_libraries(MaterialLib
    BaseLib
)
