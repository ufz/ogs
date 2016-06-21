# Source files
GET_SOURCE_FILES(SOURCES)
GET_SOURCE_FILES(SOURCES_ADSORPTION Adsorption)
set(SOURCES ${SOURCES} ${SOURCES_ADSORPTION})

add_library(MaterialsLib ${SOURCES} )

if(MSVC)
    include(GenerateExportHeader)
    generate_export_header(MaterialsLib)
endif()

target_link_libraries(MaterialsLib
    BaseLib
)
