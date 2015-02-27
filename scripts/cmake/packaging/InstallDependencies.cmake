macro(InstallDependencies TARGET INSTALL_COMPONENT)

	if(MSVC)
		set(TARGET_EXE ${EXECUTABLE_OUTPUT_PATH}/Release/${TARGET}.exe)
		set(EXCLUDE_SYSTEM 0)
	else()
		set(TARGET_EXE ${EXECUTABLE_OUTPUT_PATH}/${TARGET})
		set(EXCLUDE_SYSTEM 1)
	endif()

	if(EXISTS ${TARGET_EXE})
		include(GetPrerequisites)
		get_prerequisites(${TARGET_EXE} TARGET_DEPENDENCIES ${EXCLUDE_SYSTEM} 0 "" "")
		message(STATUS "${TARGET_EXE} dependencies:")
		foreach(DEPENDENCY ${TARGET_DEPENDENCIES})
			gp_resolve_item("/" "${DEPENDENCY}" ${TARGET_EXE}
				"/usr/local/lib;/;${VTK_DIR};/usr/lib64;" DEPENDENCY_PATH)
			get_filename_component(RESOLVED_DEPENDENCY_PATH "${DEPENDENCY_PATH}" REALPATH)
			string(TOLOWER ${DEPENDENCY} DEPENDENCY_LOWER)
				set(DEPENDENCY_PATHS ${DEPENDENCY_PATHS} ${RESOLVED_DEPENDENCY_PATH})
				message("    ${RESOLVED_DEPENDENCY_PATH}")
		endforeach()
		install(FILES ${DEPENDENCY_PATHS} DESTINATION bin COMPONENT ${INSTALL_COMPONENT})
		add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ;)
	else()
		# Run CMake after target was built to run GetPrerequisites on executable
		add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND}
			ARGS ${CMAKE_SOURCE_DIR} WORKING_DIRECTORY ${CMAKE_BINARY_DIR} VERBATIM)
	endif()

endmacro()
