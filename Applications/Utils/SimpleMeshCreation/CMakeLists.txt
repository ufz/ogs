
if(OGS_BUILD_GUI)
    add_executable(createMeshElemPropertiesFromASCRaster createMeshElemPropertiesFromASCRaster.cpp)
    set_target_properties(createMeshElemPropertiesFromASCRaster PROPERTIES FOLDER Utilities)
    target_link_libraries(createMeshElemPropertiesFromASCRaster
        ApplicationsFileIO
        MeshLib
    )
endif() # OGS_BUILD_GUI (VtkVis-target is existing)


add_executable(generateStructuredMesh generateStructuredMesh.cpp)
target_link_libraries(generateStructuredMesh MeshLib)
set_target_properties(generateStructuredMesh PROPERTIES FOLDER Utilities)

####################
### Installation ###
####################
install(TARGETS
    generateStructuredMesh
    RUNTIME DESTINATION bin
    COMPONENT Utilities
)
if(TARGET VtkVis)
    install(TARGETS
        createMeshElemPropertiesFromASCRaster
        RUNTIME DESTINATION bin
        COMPONENT Utilities
    )
endif()
