SET(CPACK_GENERATOR NSIS ZIP)
SET(CPACK_NSIS_MUI_ICON ${CMAKE_SOURCE_DIR}/scripts/packaging/ogs-de-icon.ico)
FILE(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.bmp" BACKGROUND_IMAGE)
SET(CPACK_PACKAGE_ICON ${BACKGROUND_IMAGE})
SET(CPACK_NSIS_INSTALLED_ICON_NAME ${CPACK_NSIS_MUI_ICON})
SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
SET(CPACK_NSIS_CONTACT "info@opengeosys.org")
SET(CPACK_NSIS_MODIFY_PATH OFF)
SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
SET(CPACK_NSIS_HELP_LINK "http://docs.opengeosys.org/assets/files/Documentation/User_Manual.pdf")
SET(CPACK_NSIS_MENU_LINKS
	"bin" "Executables folder"
	"http://www.opengeosys.org" "Website"
	"https://github.com/ufz/ogs" "Source code on GitHub"
	"docs" "Documentation folder"
)
