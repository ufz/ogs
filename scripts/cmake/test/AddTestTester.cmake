# Run vtk.js converter
if(VIS_FILES AND VTKJS_CONVERTER)
    execute_process(COMMAND cmake -E make_directory ${VTKJS_OUTPUT_PATH})
    foreach(FILE ${VIS_FILES})
        execute_process(
            COMMAND ${VTKJS_CONVERTER} -e -i ${BINARY_PATH}/${FILE} -o ${VTKJS_OUTPUT_PATH}
        )
    endforeach()
endif()

message(STATUS "running command checking test results: cd ${case_path} && ${TESTER_COMMAND}")

if(WIN32)
    set(TERMINAL_CMD cmd /C)
else()
    set(TERMINAL_CMD bash -c)
endif()
foreach(CMD ${TESTER_COMMAND})
    set(COMBINED_COMMAND ${COMBINED_COMMAND} COMMAND ${TERMINAL_CMD} ${CMD})
endforeach()

execute_process(
    ${COMBINED_COMMAND}
    WORKING_DIRECTORY ${case_path}
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_VARIABLE OUTPUT
)

if(NOT EXIT_CODE STREQUAL "0")
    message(FATAL_ERROR "Error exit code: ${EXIT_CODE}")
endif()
