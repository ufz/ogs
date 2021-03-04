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
    Mesh2Raster)

if(OGS_USE_NETCDF)
    list(APPEND TOOLS NetCdfConverter)
endif()

if(OGS_BUILD_GUI)
    if(Shapelib_FOUND)
        list(APPEND TOOLS ConvertSHPToGLI)
    endif()
    list(APPEND TOOLS FEFLOW2OGS)
endif()

if(Shapelib_FOUND)
    list(APPEND TOOLS Mesh2Shape)
endif()

foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} ApplicationsFileIO GitInfoLib MeshLib
        MeshGeoToolsLib tclap)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)

if(TARGET ConvertSHPToGLI)
    target_link_libraries(ConvertSHPToGLI GeoLib Qt5::Xml ${Shapelib_LIBRARIES})
endif()

if(TARGET Mesh2Shape)
    target_link_libraries(Mesh2Shape ${Shapelib_LIBRARIES})
endif()

if(TARGET NetCdfConverter)
    target_link_libraries(NetCdfConverter
                          ${NETCDF_LIBRARIES_CXX}
                          ${NETCDF_LIBRARIES_C}
                          ${HDF5_HL_LIBRARIES}
                          ${HDF5_LIBRARIES}
    )
    if(OGS_USE_CONAN AND MSVC)
        target_link_libraries(NetCdfConverter ${CONAN_LIBS})
    endif()
endif()
