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
)

if(OGS_USE_NETCDF)
    list(APPEND TOOLS NetCdfConverter)
endif()

if(OGS_BUILD_GUI)
    list(APPEND TOOLS ConvertSHPToGLI FEFLOW2OGS Mesh2Shape)
endif()

foreach(tool ${TOOLS})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(
        ${tool} ApplicationsFileIO GitInfoLib MeshLib MeshGeoToolsLib tclap
    )
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)

if(TARGET ConvertSHPToGLI)
    target_link_libraries(ConvertSHPToGLI GeoLib Qt5::Xml)
endif()

if(TARGET NetCdfConverter)
    target_link_libraries(
        NetCdfConverter ${NETCDF_LIBRARIES_CXX} ${NETCDF_LIBRARIES_C}
        ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES}
    )
    if(OGS_USE_CONAN AND MSVC)
        target_link_libraries(NetCdfConverter ${CONAN_LIBS})
    endif()
endif()
