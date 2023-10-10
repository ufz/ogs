foreach(file ${FILES_TO_DELETE})
    file(REMOVE ${file})
endforeach()

execute_process(
    COMMAND ${WRAPPER_COMMAND} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE LOG
    ERROR_VARIABLE LOG COMMAND_ECHO STDOUT
)

if(EXIT_CODE STREQUAL "0")
    if(NOT DEFINED ENV{CI} AND NOT CAT_LOG)
        file(WRITE ${LOG_FILE} "${LOG}")
    endif()
else()
    if(CAT_LOG)
        message(FATAL_ERROR "Exit code: ${EXIT_CODE}; log:\n${LOG}")
    else()
        file(WRITE ${LOG_FILE} "${LOG}")
        message(FATAL_ERROR "Exit code: ${EXIT_CODE}; log file: ${LOG_FILE}")
    endif()
endif()
