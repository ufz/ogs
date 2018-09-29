# Run vtk.js converter
if(VIS_FILES AND VTKJS_CONVERTER)
    execute_process(COMMAND cmake -E make_directory ${VTKJS_OUTPUT_PATH})
    foreach(FILE ${VIS_FILES})
        execute_process(
            COMMAND ${VTKJS_CONVERTER} -e -i ${BINARY_PATH}/${FILE} -o ${VTKJS_OUTPUT_PATH}
        )
    endforeach()
endif()

message(STATUS "running tester (glob mode: ${GLOB_MODE}): ${TESTER_COMMAND}")

if(WIN32)
    set(TERMINAL_CMD cmd /C)
else()
    set(TERMINAL_CMD bash -c)
endif()
foreach(CMD ${TESTER_COMMAND})
    if(GLOB_MODE)
        separate_arguments(CMD)
        list(GET CMD 0 GLOB)
        list(GET CMD 1 NAME_A)
        list(GET CMD 2 NAME_B)
        list(GET CMD 3 ABS_TOL)
        list(GET CMD 4 REL_TOL)
        file(GLOB FILES RELATIVE ${case_path} ${GLOB})
        list(LENGTH FILES length)
        message(STATUS "Glob expression '${GLOB}' (${NAME_A}) found ${length} files.")
        foreach(FILE ${FILES})
            execute_process(
                COMMAND ${SELECTED_DIFF_TOOL_PATH} ${case_path}/${FILE} ${BINARY_PATH}/${FILE} -a ${NAME_A} -b ${NAME_B} --abs ${ABS_TOL} --rel ${REL_TOL}
                WORKING_DIRECTORY ${case_path}
                RESULT_VARIABLE EXIT_CODE
                OUTPUT_VARIABLE OUTPUT
            )

            if(NOT EXIT_CODE STREQUAL "0")
                message(FATAL_ERROR "Error exit code: ${EXIT_CODE}\n${OUTPUT}")
            endif()
        endforeach()
    else()
        execute_process(
            COMMAND ${TERMINAL_CMD} "${CMD}"
            WORKING_DIRECTORY ${case_path}
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_VARIABLE OUTPUT
        )
        if(NOT EXIT_CODE STREQUAL "0")
            message(FATAL_ERROR "Error exit code: ${EXIT_CODE}")
        endif()
    endif()
endforeach()


