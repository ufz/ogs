GET_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES Adsorption)
APPEND_SOURCE_FILES(SOURCES SolidModels)

add_library(MaterialLib ${SOURCES})
target_link_libraries(MaterialLib
    BaseLib
)
