# Run vtk.js converter
if(VIS_FILES AND VTKJS_CONVERTER)
    execute_process(COMMAND cmake -E make_directory ${VTKJS_OUTPUT_PATH})
    foreach(file ${VIS_FILES})
        execute_process(
            COMMAND ${VTKJS_CONVERTER} -e -i ${BINARY_PATH}/${file} -o ${VTKJS_OUTPUT_PATH}
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
        message(STATUS "Glob expression '${GLOB}' (${NAME_A}) found ${LENGTH} files.")
        if(${LENGTH} EQUAL 0)
            message(FATAL_ERROR "DIFF_DATA glob expression '${GLOB}' "
                "did not match any files!")
        endif()
        foreach(file ${FILES})
            if("$ENV{HOSTNAME}" MATCHES "frontend.*")
                string(REPLACE "gpfs1" "../.." file ${file})
            endif()
            if("$ENV{HOSTNAME}" MATCHES "frontend.*")
                string(REPLACE "gpfs0" "../.." file ${file})
            endif()
            execute_process(
                COMMAND ${SELECTED_DIFF_TOOL_PATH} ${SOURCE_PATH}/${file} ${BINARY_PATH}/${file} -a ${NAME_A} -b ${NAME_B} --abs ${ABS_TOL} --rel ${REL_TOL}
                WORKING_DIRECTORY ${SOURCE_PATH}
                RESULT_VARIABLE EXIT_CODE
                OUTPUT_VARIABLE OUTPUT
                ERROR_VARIABLE OUTPUT
            )

            if(NOT EXIT_CODE STREQUAL "0")
                message(WARNING "Error exit code: ${EXIT_CODE}\n${OUTPUT}")
                set(TEST_FAILED TRUE)
            endif()
        endforeach()
    else()
        execute_process(
            COMMAND ${TERMINAL_CMD} "${cmd}"
            WORKING_DIRECTORY ${SOURCE_PATH}
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_VARIABLE OUTPUT
            ERROR_VARIABLE OUTPUT
        )
        if(NOT EXIT_CODE STREQUAL "0")
            message(WARNING "Error exit code: ${EXIT_CODE}${OUTPUT}")
            set(TEST_FAILED TRUE)
        endif()
    endif()
endforeach()
if (TEST_FAILED)
    message(FATAL_ERROR "One of the tests failed.")
endif()


