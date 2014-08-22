SET(ARGS --boundary_condition ${case_path}/${case_name}.cnd -g ${case_path}/${case_name}.gml -m ${case_path}/${case_name}.vtu)

STRING(REPLACE " " ";" WRAPPER_COMMAND ${WRAPPER_COMMAND})
EXECUTE_PROCESS(
	COMMAND ${WRAPPER_COMMAND} ${executable} ${ARGS}
	WORKING_DIRECTORY ${case_path}
	RESULT_VARIABLE EXIT_CODE
)

IF(NOT EXIT_CODE STREQUAL "0")
	MESSAGE(FATAL_ERROR "Test wrapper exited with code: ${EXIT_CODE}")
ENDIF()
