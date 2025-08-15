# IMPORTANT: multiple arguments in a single variable have to be in list notation (;)
# and have to be quoted when passed
# "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
execute_process(
    COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${UV_RUN_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE LOG
    ERROR_VARIABLE LOG
    ECHO_OUTPUT_VARIABLE
    ECHO_ERROR_VARIABLE
    COMMAND_ECHO STDOUT
)

set(SAVE_LOG true)
set(TEST_LOG_DIR "${LOG_ROOT}")

if (TEST_COMMAND_IS_EXPECTED_TO_SUCCEED)
    if (EXIT_CODE STREQUAL "0")
        # expected: success, actual: success
        if (DEFINED ENV{CI} AND NOT OGS_CI_ALWAYS_SAVE_LOG_FILE_TO_ARTIFACTS)
            set(SAVE_LOG false)
        endif()
    else()
        # expected: success, actual: failure
        # use default settings
    endif()
else()
    if (EXIT_CODE STREQUAL "0")
        # expected: failure, actual: success
        if (DEFINED ENV{CI})
            set(TEST_LOG_DIR "${LOG_ROOT}/command_succeeded_but_was_expected_to_fail")
        endif()
    else()
        # expected: failure, actual: failure
        if (DEFINED ENV{CI})
            set(TEST_LOG_DIR "${LOG_ROOT}/command_failed_as_expected")
        endif()
    endif()
endif()

set(LOG_FILE "${TEST_LOG_DIR}/${LOG_FILE_BASENAME}")

if (SAVE_LOG)
    if(NOT EXISTS "${TEST_LOG_DIR}")
        file(MAKE_DIRECTORY "${TEST_LOG_DIR}")
    endif()

    file(WRITE "${LOG_FILE}" "${LOG}")
endif()

if(NOT EXIT_CODE STREQUAL "0")
    message(FATAL_ERROR "Exit code: ${EXIT_CODE}; log file: ${LOG_FILE}")
endif()
