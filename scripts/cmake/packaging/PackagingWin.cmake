set(CMAKE_INSTALL_UCRT_LIBRARIES ON)
set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
if(NOT OGS_INSTALL_DEPENDENCIES)
    include(InstallRequiredSystemLibraries)
endif()
set(CPACK_GENERATOR ZIP)
set(CPACK_NSIS_MUI_ICON ${PROJECT_SOURCE_DIR}/scripts/packaging/ogs-de-icon.ico)
set(CPACK_PACKAGE_ICON ${PROJECT_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.bmp)
set(CPACK_NSIS_INSTALLED_ICON_NAME ${CPACK_NSIS_MUI_ICON})
set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
set(CPACK_NSIS_CONTACT "info@opengeosys.org")
set(CPACK_NSIS_MODIFY_PATH OFF)
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
set(CPACK_NSIS_HELP_LINK "http://docs.opengeosys.org/assets/files/Documentation/User_Manual.pdf")
set(CPACK_NSIS_MENU_LINKS
    "bin" "Executables folder"
    "http://www.opengeosys.org" "Website"
    "https://github.com/ufz/ogs" "Source code on GitHub"
)

if(OGS_USE_CONAN)
    file(GLOB MATCHED_FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*.dll")
    install(FILES ${MATCHED_FILES} DESTINATION bin)
endif()

if(OGS_BUILD_GUI)
    install_qt5_plugin("Qt5::QWindowsIntegrationPlugin" QT_PLUGINS)
    file(WRITE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/qt.conf"
        "[Paths]\nPlugins = ../${_qt_plugin_dir}\n")
    install(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/qt.conf"
        DESTINATION bin)
endif()
