if(OGS_USE_MPI)
    return()
endif()

set(TOOLS generateStructuredMesh)
if(OGS_BUILD_GUI)
    list(APPEND TOOLS createMeshElemPropertiesFromASCRaster)
endif()

foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} ApplicationsFileIO GitInfoLib MeshLib tclap)
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
