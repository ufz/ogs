INCLUDE(packaging/PackagingMacros)

#### Packaging setup ####
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OGS-6 THM/C Simulator")
SET(CPACK_PACKAGE_VENDOR "OpenGeoSys Community (http://www.opengeosys.org)")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OGS-${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
# SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_PACKAGE_VERSION_MAJOR "${OGS_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${OGS_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${OGS_VERSION_PATCH}")
IF(APPLE)
	SET(CPACK_PACKAGE_FILE_NAME "ogs-${OGS_VERSION}-OSX-x${BITS}")
	SET(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
ELSE()
	SET(CPACK_PACKAGE_FILE_NAME "ogs-${OGS_VERSION}-${CMAKE_SYSTEM}-x${BITS}")
ENDIF()

#SET(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE 1)

IF (WIN32)
	INCLUDE (packaging/PackagingWin)
ENDIF()
IF(UNIX)
	INCLUDE (packaging/PackagingLinux)
ENDIF()
IF(APPLE)
	INCLUDE (packaging/PackagingMac)
ENDIF()

# Download additional content
IF(OGS_DOWNLOAD_ADDITIONAL_CONTENT)
	DownloadAdditionalFilesForPackaging(
		URLS http://docs.opengeosys.org/assets/releases/head/docs/DataExplorer-Manual.pdf
		     http://docs.opengeosys.org/assets/releases/head/docs/User_Manual.pdf
		     http://docs.opengeosys.org/assets/releases/head/docs/Theory_Manual.pdf
		DESTINATION docs
		PACKAGE_GROUP ogs_docs
	)

	IF(WIN32)
		DownloadAdditionalFilesForPackaging(
			URLS http://docs.opengeosys.org/assets/releases/head/win/OGSFileConverter.exe
			DESTINATION bin
			EXECUTABLE TRUE
			PACKAGE_GROUP ogs_converter
		)
	ENDIF()
	IF(APPLE)
		DownloadAdditionalFilesForPackaging(
			URLS http://docs.opengeosys.org/assets/releases/head/mac/OGSFileConverter
			DESTINATION bin
			EXECUTABLE TRUE
			PACKAGE_GROUP ogs_converter
		)
	ENDIF()
ENDIF()

INCLUDE (CPack)

cpack_add_component_group(Applications
	DISPLAY_NAME Applications
	DESCRIPTION "OpenGeoSys applications"
	EXPANDED
	BOLD_TITLE
)

cpack_add_component_group(Utilities
	DISPLAY_NAME Utilities
	DESCRIPTION "OpenGeoSys utilities"
	EXPANDED
)

cpack_add_component(ogs_extras
	DISPLAY_NAME "Extra tools"
	DESCRIPTION "Miscellaneous tools."
	GROUP Utilities
)

cpack_add_component(ogs_docs
	DISPLAY_NAME "Documentation"
	DESCRIPTION "PDF documentation."
	GROUP Utilities
)
