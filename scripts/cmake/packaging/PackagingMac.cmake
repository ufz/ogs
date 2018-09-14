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

SET(CMAKE_INSTALL_RPATH "@executable_path;@executable_path/../${CMAKE_INSTALL_LIBDIR}")

if(OGS_BUILD_GUI)
    install_qt5_plugin("Qt5::QCocoaIntegrationPlugin" QT_PLUGINS)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        "[Paths]\nPlugins = ../${_qt_plugin_dir}\n")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        DESTINATION bin COMPONENT ogs_gui)
endif()
