if(APPLE)
	find_program(REALPATH_TOOL_PATH grealpath)
elseif(UNIX)
	find_program(REALPATH_TOOL_PATH realpath)
else()
	return()
endif()
if(NOT REALPATH_TOOL_PATH)
	return()
endif()

add_custom_target(archive-data
	bash ${CMAKE_SOURCE_DIR}/scripts/packaging/archive-testdata.sh
	DEPENDS data
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Packaging testdata to ogs6-data.tar.gz" VERBATIM
)
