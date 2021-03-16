# IMPORTANT: multiple arguments in one variables have to be in list notation (;)
# and have to be quoted when passed
# "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
foreach(file ${FILES_TO_DELETE})
    file(REMOVE ${BINARY_PATH}/${file})
endforeach()

string(REPLACE ";" " " CMD_STRING
               "cd ${WORKING_DIRECTORY} && ${WRAPPER_COMMAND} "
               "${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}"
)
message(STATUS "Test command cleaned:\n${CMD_STRING}")

execute_process(
    COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${WORKING_DIRECTORY}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE OUTPUT
    ERROR_VARIABLE OUTPUT
)

if(NOT EXIT_CODE STREQUAL "0")
    message(
        FATAL_ERROR "Test wrapper exited with code: ${EXIT_CODE}\n${OUTPUT}"
    )
endif()
