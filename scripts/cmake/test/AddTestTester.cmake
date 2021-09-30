set(_exec_process_args "")
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
    set(_exec_process_args ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE)
endif()

# Run vtk.js converter
if(VIS_FILES AND VTKJS_CONVERTER)
    execute_process(COMMAND cmake -E make_directory ${VTKJS_OUTPUT_PATH})
    foreach(file ${VIS_FILES})
        execute_process(
            COMMAND ${VTKJS_CONVERTER} -e -i ${BINARY_PATH}/${file} -o
                    ${VTKJS_OUTPUT_PATH}
        )
    endforeach()
endif()

message(STATUS "running tester (glob mode: ${GLOB_MODE}): ${TESTER_COMMAND}")

if(WIN32)
    set(TERMINAL_CMD cmd /C)
else()
    set(TERMINAL_CMD bash -c)
endif()
set(TEST_FAILED FALSE)
set(_counter 0)
foreach(cmd ${TESTER_COMMAND})
    if(GLOB_MODE)
        # cmake-lint: disable=E1120
        separate_arguments(cmd)
        list(GET cmd 0 GLOB)
        list(GET cmd 1 NAME_A)
        list(GET cmd 2 NAME_B)
        list(GET cmd 3 ABS_TOL)
        list(GET cmd 4 REL_TOL)
        file(GLOB FILES RELATIVE ${SOURCE_PATH} ${GLOB})
        list(LENGTH FILES LENGTH)
        message(
            STATUS
                "Glob expression '${GLOB}' (${NAME_A}) found ${LENGTH} files."
        )
        if(${LENGTH} EQUAL 0)
            message(FATAL_ERROR "DIFF_DATA glob expression '${GLOB}' "
                                "did not match any files!"
            )
        endif()
        foreach(file ${FILES})
            math(EXPR _counter "${_counter}+1")
            set(LOG_FILE ${LOG_FILE_BASE}-${_counter}.txt)
            if("$ENV{HOSTNAME}" MATCHES "frontend.*")
                string(REPLACE "gpfs1" "../.." file ${file})
            endif()
            if("$ENV{HOSTNAME}" MATCHES "frontend.*")
                string(REPLACE "gpfs0" "../.." file ${file})
            endif()
            execute_process(
                COMMAND
                    ${SELECTED_DIFF_TOOL_PATH} ${SOURCE_PATH}/${file}
                    ${BINARY_PATH}/${file} -a ${NAME_A} -b ${NAME_B} --abs
                    ${ABS_TOL} --rel ${REL_TOL}
                WORKING_DIRECTORY ${SOURCE_PATH}
                RESULT_VARIABLE EXIT_CODE
                OUTPUT_VARIABLE OUTPUT
                ERROR_VARIABLE OUTPUT ${_exec_process_args}
            )
            if(NOT EXIT_CODE STREQUAL "0")
                file(WRITE ${LOG_FILE} ${OUTPUT})
                message(
                    WARNING "Exit code: ${EXIT_CODE}; log file: ${LOG_FILE}"
                )
                set(TEST_FAILED TRUE)
            endif()
        endforeach()
    else()
        math(EXPR _counter "${_counter}+1")
        set(LOG_FILE ${LOG_FILE_BASE}-${_counter}.txt)
        execute_process(
            COMMAND ${TERMINAL_CMD} "${cmd}"
            WORKING_DIRECTORY ${SOURCE_PATH}
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_VARIABLE OUTPUT
            ERROR_VARIABLE OUTPUT ${_exec_process_args}
        )
        if(NOT EXIT_CODE STREQUAL "0")
            file(WRITE ${LOG_FILE} ${OUTPUT})
            message(WARNING "Exit code: ${EXIT_CODE}; log file: ${LOG_FILE}")
            set(TEST_FAILED TRUE)
        endif()
    endif()
endforeach()
if(TEST_FAILED)
    message(FATAL_ERROR "One of the tests failed.")
endif()
