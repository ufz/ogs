# IMPORTANT: multiple arguments in one variables have to be in list notation (;)
# and have to be quoted when passed "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
execute_process(
	COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
	WORKING_DIRECTORY ${case_path}
	RESULT_VARIABLE EXIT_CODE
)

if(NOT EXIT_CODE STREQUAL "0")
	message(FATAL_ERROR "Test wrapper exited with code: ${EXIT_CODE}")
endif()
