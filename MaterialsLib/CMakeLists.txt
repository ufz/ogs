# Source files
GET_SOURCE_FILES(SOURCES_MATERIALSLIB)

add_subdirectory( Adsorption )

add_library(MaterialsLib ${SOURCES_MATERIALSLIB} )
