# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)
append_source_files(SOURCES SolidModels)

append_source_files(SOURCES Fluid/Density)
append_source_files(SOURCES Solid/Density)

add_library(MaterialLib ${SOURCES} )
target_link_libraries(MaterialLib
    BaseLib
)
