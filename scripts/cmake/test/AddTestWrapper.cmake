# IMPORTANT: multiple arguments in one variables have to be in list notation (;)
# and have to be quoted when passed "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
foreach(FILE ${FILES_TO_DELETE})
    file(REMOVE ${BINARY_PATH}/${FILE})
endforeach()

# taken from https://stackoverflow.com/a/7216542
function(JOIN VALUES GLUE OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

JOIN("${WRAPPER_ARGS}" " " WRAPPER_ARGS_STR)
JOIN("${EXECUTABLE_ARGS}" " " EXECUTABLE_ARGS_STR)

message(STATUS "running command generating test results: cd ${case_path} && ${WRAPPER_COMMAND} ${WRAPPER_ARGS_STR} ${EXECUTABLE} ${EXECUTABLE_ARGS_STR}")
execute_process(
    COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${case_path}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_FILE ${STDOUT_FILE_PATH}
)

if(NOT EXIT_CODE STREQUAL "0")
    message(FATAL_ERROR "Test wrapper exited with code: ${EXIT_CODE}")
endif()
