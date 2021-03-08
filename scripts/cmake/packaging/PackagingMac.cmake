set(CPACK_GENERATOR TGZ)

option(OGS_PACK_DMG "package targets creates a .dmg disk image instead of .tar.gz" FALSE)
if(OGS_PACK_DMG)
    if(NOT OGS_BUILD_GUI)
        message(WARNING "OGS_PACK_DMG requires OGS_BUILD_GUI=ON!")
    endif()
    set(CPACK_GENERATOR DragNDrop)
endif()
set(CPACK_DMG_FORMAT "UDBZ")

# See http://stackoverflow.com/a/16662169/80480 how to create the DS_Store file.
set(CPACK_DMG_BACKGROUND_IMAGE ${PROJECT_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.png)
set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/scripts/packaging/.DS_Store)

if(OGS_USE_CONAN)
    file(GLOB MATCHED_FILES LIST_DIRECTORIES false "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/*.dylib*")
    install(FILES ${MATCHED_FILES} DESTINATION lib)

    # macOS frameworks are directories, exclude header files
    file(GLOB MATCHED_DIRECTORIES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*.framework")
    install(DIRECTORY ${MATCHED_DIRECTORIES} DESTINATION bin
        PATTERN "Headers" EXCLUDE)
endif()

if(OGS_BUILD_GUI)
    install_qt5_plugin("Qt5::QCocoaIntegrationPlugin" QT_PLUGINS)
    file(WRITE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/qt.conf"
        "[Paths]\nPlugins = ../${_qt_plugin_dir}\n")
    install(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/qt.conf"
        DESTINATION bin)
endif()
