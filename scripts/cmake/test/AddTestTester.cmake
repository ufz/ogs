#message("tester: ${TESTER_COMMAND}")
#string(REPLACE " " ";" TESTER_COMMAND ${TESTER_COMMAND})
#set(list ${TESTER_COMMAND})
#message("tester: ${list}")


execute_process(
	COMMAND bash -c ${TESTER_COMMAND}
	WORKING_DIRECTORY ${case_path}
	RESULT_VARIABLE EXIT_CODE
)

if(NOT EXIT_CODE STREQUAL "0")
	message(FATAL_ERROR "Error exit code: ${EXIT_CODE}")
endif()
