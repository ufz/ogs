# This set MSVC_INSTALL_PATHS cache variable to a list of possible Visual Studio
# install directories.
# Usage:
#
#  INCLUDE(MSVCPaths)
#    FIND_PROGRAM(DUMPBIN_TOOL_PATH dumpbin DOC "Windows dependency analysis tool"
#      PATHS MSVC_INSTALL_PATHS PATH_SUFFIXES VC/bin)

IF(MSVC)
	IF(MSVC_VERSION EQUAL 1700)
		SET(MSVC_NUMBER 11.0)
	ELSEIF(MSVC_VERSION EQUAL 1800)
		SET(MSVC_NUMBER 12.0)
	ENDIF()
	GET_FILENAME_COMPONENT(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\${MSVC_NUMBER}\\Setup\\VS;ProductDir]" REALPATH)
	GET_FILENAME_COMPONENT(VS_EXPRESS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VCExpress\\${MSVC_NUMBER}\\Setup\\VS;ProductDir]" REALPATH)
	SET(MSVC_INSTALL_PATHS
		${VS_DIR} ${VS_EXPRESS_DIR}
		"$ENV{ProgramFiles}/Microsoft Visual Studio ${MSVC_NUMBER}"
		"$ENV{ProgramFiles(x86)}/Microsoft Visual Studio ${MSVC_NUMBER}"
		CACHE STRING "" FORCE)
ENDIF()
