# IMPORTANT: multiple arguments in one variables have to be in list notation (;)
# and have to be quoted when passed
# "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
execute_process(
    COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE LOG
    ERROR_VARIABLE LOG
    ECHO_OUTPUT_VARIABLE
    ECHO_ERROR_VARIABLE
    COMMAND_ECHO STDOUT
)

if(EXIT_CODE STREQUAL "0" AND NOT DEFINED ENV{CI})
    file(WRITE ${LOG_FILE} "${LOG}")
elseif(NOT EXIT_CODE STREQUAL "0")
    file(WRITE ${LOG_FILE} "${LOG}")
    message(FATAL_ERROR "Exit code: ${EXIT_CODE}; log file: ${LOG_FILE}")
endif()
