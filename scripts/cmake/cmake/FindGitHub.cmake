# - Find GitHub for Windows
#
#   GITHUB_FOUND    - Was GitHub for Windows found
#   GITHUB_BIN_DIR  - Path to the bin-directory where useful bash tools can be found
#
# Example usage:
#   FIND_PACKAGE(GitHub)
#   FIND_PROGRAM(BASH_TOOL_PATH bash HINTS ${GITHUB_BIN_DIR} DOC "The bash executable")

IF(WIN32 AND NOT GITHUB_FOUND)

	# Check install Path
	FIND_PATH(
		GITHUB_DIR
		shell.ps1
		PATHS $ENV{LOCALAPPDATA}/GitHub
		NO_DEFAULT_PATH
	)

	IF(GITHUB_DIR)

		FILE(TO_NATIVE_PATH ${GITHUB_DIR} GITHUB_WIN_DIR)
		EXECUTE_PROCESS (
			COMMAND cmd /c "cd ${GITHUB_WIN_DIR}/PortableGit* & cd"
			OUTPUT_VARIABLE PORTABLE_GIT_WIN_DIR
		)

		IF(PORTABLE_GIT_WIN_DIR)
			STRING(STRIP ${PORTABLE_GIT_WIN_DIR} PORTABLE_GIT_WIN_DIR)
			FILE(TO_CMAKE_PATH ${PORTABLE_GIT_WIN_DIR} PORTABLE_GIT_DIR)
			SET(GITHUB_FOUND ON CACHE BOOL "Was GitHub for Windows found?")
			SET(GITHUB_BIN_DIR ${PORTABLE_GIT_DIR}/bin CACHE PATH "The path to the GitHub for Windows binaries.")
			MESSAGE(STATUS "GitHub for Windows found.")
		ENDIF()

	ENDIF()
ENDIF()