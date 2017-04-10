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
    target_link_libraries(${TOOL} ApplicationsFileIO GitInfoLib MeshLib)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin COMPONENT Utilities)

if(TARGET ConvertSHPToGLI)
    target_link_libraries(ConvertSHPToGLI GeoLib Qt5::Xml ${Shapelib_LIBRARIES})
endif()

add_executable(OGS2VTK OGS2VTK.cpp)
set_target_properties(OGS2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(OGS2VTK MeshLib)
ADD_VTK_DEPENDENCY(OGS2VTK)

add_executable(VTK2OGS VTK2OGS.cpp)
set_target_properties(VTK2OGS PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2OGS MeshLib)
ADD_VTK_DEPENDENCY(VTK2OGS)

add_executable(VTK2TIN VTK2TIN.cpp)
set_target_properties(VTK2TIN PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2TIN MeshLib)
ADD_VTK_DEPENDENCY(VTK2TIN)

add_executable(TIN2VTK TIN2VTK.cpp)
set_target_properties(TIN2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(TIN2VTK MeshLib)
ADD_VTK_DEPENDENCY(TIN2VTK)

add_executable(NetCdfConverter NetCdfConverter.cpp)
set_target_properties(NetCdfConverter PROPERTIES FOLDER Utilities)
target_link_libraries(NetCdfConverter MeshLib)
ADD_VTK_DEPENDENCY(NetCdfConverter)

####################
### Installation ###
####################
install(TARGETS generateMatPropsFromMatID GMSH2OGS OGS2VTK VTK2OGS VTK2TIN
    RUNTIME DESTINATION bin COMPONENT ogs_converter)

if(Qt5XmlPatterns_FOUND)
    if(Shapelib_FOUND)
        install(TARGETS ConvertSHPToGLI
            RUNTIME DESTINATION bin COMPONENT ogs_converter)
    endif()
    install(TARGETS  FEFLOW2OGS convertGEO
        RUNTIME DESTINATION bin COMPONENT ogs_converter)
endif()

