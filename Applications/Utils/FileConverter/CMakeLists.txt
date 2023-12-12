set(TOOLS
    convertGEO
    generateMatPropsFromMatID
    GMSH2OGS
    OGS2VTK
    VTK2OGS
    VTK2TIN
    TIN2VTK
    TecPlotTools
    GocadSGridReader
    GocadTSurfaceReader
    Mesh2Raster
    PVD2XDMF
    Raster2ASC
)

if(OGS_USE_NETCDF)
    list(APPEND TOOLS NetCdfConverter)
endif()

if(OGS_BUILD_GUI)
    list(APPEND TOOLS ConvertSHPToGLI FEFLOW2OGS Mesh2Shape)
endif()

foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(
        ${tool} ApplicationsFileIO GitInfoLib MeshLib
                MeshToolsLib MeshGeoToolsLib tclap
    )
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)

if(TARGET ConvertSHPToGLI)
    target_link_libraries(ConvertSHPToGLI GeoLib Qt5::Xml)
endif()

if(TARGET NetCdfConverter)
    target_link_libraries(NetCdfConverter netcdf-cxx4)
endif()
