# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)
append_source_files(SOURCES SolidModels)
append_source_files(SOURCES FractureModels)

append_source_files(SOURCES Fluid)
append_source_files(SOURCES Fluid/Density)

add_library(MaterialLib ${SOURCES} )
target_link_libraries(MaterialLib
    BaseLib
)
