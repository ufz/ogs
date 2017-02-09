if(WIN32)
    execute_process(
        COMMAND cmd /C ${TESTER_COMMAND}
        WORKING_DIRECTORY ${case_path}
        RESULT_VARIABLE EXIT_CODE
        OUTPUT_VARIABLE OUTPUT
  )
else()
    execute_process(
        COMMAND bash -c ${TESTER_COMMAND}
        WORKING_DIRECTORY ${case_path}
        RESULT_VARIABLE EXIT_CODE
        OUTPUT_VARIABLE OUTPUT
    )
endif()

if(NOT EXIT_CODE STREQUAL "0")
    message(FATAL_ERROR "Error exit code: ${EXIT_CODE}")
endif()

# Run vtk.js converter
if(VIS_FILES)
    execute_process(COMMAND cmake -E make_directory ${VTKJS_OUTPUT_PATH})
endif()
foreach(FILE ${VIS_FILES})
    execute_process(
        COMMAND ${VTKJS_CONVERTER} -e -i ${BINARY_PATH}/${FILE} -o ${VTKJS_OUTPUT_PATH}
    )
endforeach()
