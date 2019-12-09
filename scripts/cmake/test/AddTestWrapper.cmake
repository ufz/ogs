# IMPORTANT: multiple arguments in one variables have to be in list notation (;)
# and have to be quoted when passed "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
foreach(FILE ${FILES_TO_DELETE})
    file(REMOVE ${BINARY_PATH}/${FILE})
endforeach()

# Create Python virtual environment and install packages
set(PIP .venv/bin/pip)
if(WIN32)
    set(PIP .venv/Scripts/pip.exe)
endif()
if(EXISTS ${SOURCE_PATH}/requirements.txt AND NOT EXISTS ${BINARY_PATH}/${PIP})
    message(STATUS "Generating Python virtual environment...")
    execute_process(
        COMMAND virtualenv .venv
        WORKING_DIRECTORY ${BINARY_PATH})
endif()
if(EXISTS ${SOURCE_PATH}/requirements.txt)
    execute_process(
        COMMAND ${PIP} install -r ${SOURCE_PATH}/requirements.txt
        WORKING_DIRECTORY ${BINARY_PATH})
endif()

execute_process(
    COMMAND ${WRAPPER_COMMAND} ${WRAPPER_ARGS} ${EXECUTABLE} ${EXECUTABLE_ARGS}
    WORKING_DIRECTORY ${SOURCE_PATH}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE OUTPUT
    ERROR_VARIABLE OUTPUT
)

if(NOT EXIT_CODE STREQUAL "0")
    message(FATAL_ERROR "Test wrapper exited with code: ${EXIT_CODE}\n${OUTPUT}")
endif()
