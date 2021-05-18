get_source_files(SOURCES)
append_source_files(SOURCES Permeability)

ogs_add_library(MaterialLib_FractureModels ${SOURCES})
target_link_libraries(
    MaterialLib_FractureModels PUBLIC MathLib PRIVATE ParameterLib NumLib
)
