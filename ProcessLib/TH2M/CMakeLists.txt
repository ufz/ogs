get_source_files(SOURCES PhaseTransitionModels)

ogs_add_library(TH2M ${SOURCES})
target_link_libraries(TH2M PUBLIC ProcessLib PRIVATE ParameterLib)
