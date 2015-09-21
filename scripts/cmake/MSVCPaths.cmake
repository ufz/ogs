# This set MSVC_INSTALL_PATHS cache variable to a list of possible Visual Studio
# install directories.
# Usage:
#
#  include(MSVCPaths)
#    find_program(DUMPBIN_TOOL_PATH dumpbin DOC "Windows dependency analysis tool"
#      PATHS MSVC_INSTALL_PATHS PATH_SUFFIXES VC/bin)

if(MSVC)
	if(MSVC_VERSION EQUAL 1700)
		set(MSVC_NUMBER 11.0)
	elseif(MSVC_VERSION EQUAL 1800)
		set(MSVC_NUMBER 12.0)
	endif()
	get_filename_component(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\${MSVC_NUMBER}\\Setup\\VS;ProductDir]" REALPATH)
	get_filename_component(VS_EXPRESS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VCExpress\\${MSVC_NUMBER}\\Setup\\VS;ProductDir]" REALPATH)

	set(X86_TMP "ProgramFiles(x86)")
	set(MSVC_INSTALL_PATHS
		${VS_DIR} ${VS_EXPRESS_DIR}
		"$ENV{ProgramFiles}/Microsoft\ Visual\ Studio\ ${MSVC_NUMBER}"
		"$ENV{${X86_TMP}}/Microsoft\ Visual\ Studio\ ${MSVC_NUMBER}"
		CACHE STRING "" FORCE)
endif()
