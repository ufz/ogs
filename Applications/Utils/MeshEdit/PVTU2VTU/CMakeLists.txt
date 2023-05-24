ogs_add_executable(
    pvtu2vtu
    PVTU2VTU.cpp
    ${PROJECT_SOURCE_DIR}/Applications/Utils/ModelPreparation/PartitionMesh/IntegrationPointDataTools.cpp
)
target_link_libraries(pvtu2vtu GitInfoLib tclap MeshToolsLib NumLib)
install(TARGETS pvtu2vtu RUNTIME DESTINATION bin)
