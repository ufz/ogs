if(OGS_BUILD_PROCESS_LIE)
    ogs_add_executable(postLIE postLIE.cpp)
    target_link_libraries(postLIE GitInfoLib LIECommon tclap)
    install(TARGETS postLIE RUNTIME DESTINATION bin)
endif()

if(OGS_BUILD_GUI)
    ogs_add_executable(Raster2PointCloud Raster2PointCloud.cpp)
    target_link_libraries(
        Raster2PointCloud ApplicationsFileIO BaseLib GitInfoLib tclap
        VtkVisFilter
    )
    install(TARGETS Raster2PointCloud RUNTIME DESTINATION bin)
endif()
