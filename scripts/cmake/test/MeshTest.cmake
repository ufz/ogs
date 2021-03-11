#
# MeshTest
# -------
#
# ~~~
# MeshTest(
#   NAME <name of the the test>
#   PATH <working directory> # relative to SourceDir/Tests/Data
#   EXECUTABLE <executable target> # optional, defaults to ogs
#   EXECUTABLE_ARGS <arguments>
#   WRAPPER <time|memcheck|callgrind|mpirun> # optional
#   WRAPPER_ARGS <arguments> # optional
#   REQUIREMENTS # optional simple boolean expression which has to be true to
#                  enable the test, e.g.
#                  OGS_USE_PETSC AND (FOO OR BAR)
#   RUNTIME <in seconds> # optional for optimizing ctest duration
#                          values should be taken from envinf1 serial job
# )
# ~~~
# cmake-lint: disable=C0103
function (MeshTest)
    if(NOT OGS_BUILD_TESTING)
        return()
    endif()
    # parse arguments
    set(options NONE)
    set(oneValueArgs EXECUTABLE PATH NAME WRAPPER RUNTIME WORKING_DIRECTORY)
    set(multiValueArgs EXECUTABLE_ARGS DATA DIFF_DATA WRAPPER_ARGS REQUIREMENTS)
    cmake_parse_arguments(MeshTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(MeshTest_SOURCE_PATH "${Data_SOURCE_DIR}/${MeshTest_PATH}")
    set(MeshTest_BINARY_PATH "${Data_BINARY_DIR}/${MeshTest_PATH}")
    file(MAKE_DIRECTORY ${MeshTest_BINARY_PATH})
    file(TO_NATIVE_PATH "${MeshTest_BINARY_PATH}" MeshTest_BINARY_PATH_NATIVE)
    set(MeshTest_STDOUT_FILE_PATH "${MeshTest_BINARY_PATH}/${MeshTest_NAME}_stdout.log")

    # set defaults
    if (NOT DEFINED MeshTest_REQUIREMENTS)
        set(MeshTest_REQUIREMENTS TRUE)
    endif()
    if (NOT DEFINED MeshTest_RUNTIME)
        set(MeshTest_RUNTIME 1)
    endif()
    if(NOT DEFINED MeshTest_WORKING_DIRECTORY)
        set(MeshTest_WORKING_DIRECTORY ${MeshTest_BINARY_PATH})
    endif()

    # --- Implement wrappers ---
    # check requirements, disable if not met
    if(${MeshTest_REQUIREMENTS})
        # message(STATUS "Enabling test ${MeshTest_NAME}.")
    else()
        set(DISABLED_TESTS_LOG "${DISABLED_TESTS_LOG}\nRequirement ${MeshTest_REQUIREMENTS} not met! Disabling test ${MeshTest_NAME}." CACHE INTERNAL "")
        return()
    endif()

    if(MeshTest_WRAPPER STREQUAL "time")
        if(TIME_TOOL_PATH)
            set(WRAPPER_COMMAND time)
        else()
            set(DISABLED_TESTS_LOG "${DISABLED_TESTS_LOG}\nDisabling time wrapper for ${MeshTest_NAME} as time exe was not found!" CACHE INTERNAL "")
            set(MeshTest_WRAPPER_ARGS "")
        endif()
    elseif(MeshTest_WRAPPER STREQUAL "memcheck")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${MeshTest_SOURCE_PATH}/${MeshTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
            set(tester memcheck)
        else()
            set(DISABLED_TESTS_LOG "${DISABLED_TESTS_LOG}\nDisabling memcheck wrapper for ${MeshTest_NAME} as memcheck exe was not found!" CACHE INTERNAL "")
            set(MeshTest_WRAPPER_ARGS "")
        endif()
    elseif(MeshTest_WRAPPER STREQUAL "callgrind")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
            unset(tester)
        else()
            set(DISABLED_TESTS_LOG "${DISABLED_TESTS_LOG}\nDisabling callgrind wrapper for ${MeshTest_NAME} as callgrind exe was not found!" CACHE INTERNAL "")
            set(MeshTest_WRAPPER_ARGS "")
        endif()
    elseif(MeshTest_WRAPPER STREQUAL "mpirun")
        if(MPIRUN_TOOL_PATH)
            set(WRAPPER_COMMAND ${MPIRUN_TOOL_PATH})
        else()
            message(STATUS "ERROR: mpirun was not found but is required for ${MeshTest_NAME}!")
            return()
        endif()
    endif()

    # --- Implement testers ---
    if(NOT MeshTest_DIFF_DATA)
        message(FATAL_ERROR "MeshTest(): ${MeshTest_NAME} - no DIFF_DATA given!")
    endif()

    string(LENGTH "${MeshTest_DIFF_DATA}" DIFF_DATA_LENGTH)
    if(${DIFF_DATA_LENGTH} GREATER 7500)
        message(FATAL_ERROR "${MeshTest_NAME}: DIFF_DATA to long! Consider using regex-syntax: TODO")
    endif()

    set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)

    list(LENGTH MeshTest_DIFF_DATA DiffDataLength)
    math(EXPR DiffDataLengthMod3 "${DiffDataLength} % 3")
    if (${DiffDataLengthMod3} EQUAL 0)
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 4)
            list(GET MeshTest_DIFF_DATA "${DiffDataIndex}" REFERENCE_VTK_FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET MeshTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET MeshTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOLERANCE)

            list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} -m \
            ${MeshTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
            ${MeshTest_BINARY_PATH}/${VTK_FILE} \
            --abs ${ABS_TOLERANCE}")
        endforeach()
    else ()
        message(FATAL_ERROR "The number of diff data arguments must be a
        multiple of three: expected.vtu output.vtu absolute_tolerance.")
    endif()

    ## -----------
    if(TARGET ${MeshTest_EXECUTABLE})
        set(MeshTest_EXECUTABLE_PARSED $<TARGET_FILE:${MeshTest_EXECUTABLE}>)
    else()
        set(MeshTest_EXECUTABLE_PARSED ${MeshTest_EXECUTABLE})
    endif()

    set(FILES_TO_DELETE "")
    list(APPEND FILES_TO_DELETE "${MeshTest_STDOUT_FILE_PATH}")
    foreach(ITEM ${MeshTest_DIFF_DATA})
        if(ITEM MATCHES "^.*\.(vtu|vtk)$")
            list(APPEND FILES_TO_DELETE "${ITEM}")
        endif()
    endforeach()

    # Run the wrapper
    if(DEFINED MeshTest_WRAPPER)
        set(MeshTest_WRAPPER_STRING "-${MeshTest_WRAPPER}")
    endif()
    set(TEST_NAME "${MeshTest_EXECUTABLE}-${MeshTest_NAME}${MeshTest_WRAPPER_STRING}")
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
        -DEXECUTABLE=${MeshTest_EXECUTABLE_PARSED}
        "-DEXECUTABLE_ARGS=${MeshTest_EXECUTABLE_ARGS}" # Quoted because passed as list
                                                        # see https://stackoverflow.com/a/33248574/80480
        -DBINARY_PATH=${MeshTest_BINARY_PATH}
        -DWRAPPER_COMMAND=${WRAPPER_COMMAND}
        "-DWRAPPER_ARGS=${MeshTest_WRAPPER_ARGS}"
        "-DFILES_TO_DELETE=${FILES_TO_DELETE}"
        -DSTDOUT_FILE_PATH=${MeshTest_STDOUT_FILE_PATH}
        -DWORKING_DIRECTORY=${MeshTest_WORKING_DIRECTORY}
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )
    set_tests_properties(${TEST_NAME} PROPERTIES COST ${MeshTest_RUNTIME})

    if(TARGET ${MeshTest_EXECUTABLE})
        add_dependencies(ctest ${MeshTest_EXECUTABLE})
        add_dependencies(ctest-large ${MeshTest_EXECUTABLE})
    endif()

    # Run the tester
    set(MeshTest_TESTER "vtkdiff")
    set(TESTER_NAME "${TEST_NAME}-${MeshTest_TESTER}")
    add_test(
        NAME ${TESTER_NAME}
        COMMAND ${CMAKE_COMMAND}
        -DSOURCE_PATH=${MeshTest_SOURCE_PATH}
        -DBINARY_PATH=${${MeshTest_BINARY_PATH}}
        -DSELECTED_DIFF_TOOL_PATH=${SELECTED_DIFF_TOOL_PATH}
        "-DTESTER_COMMAND=${TESTER_COMMAND}"
        -DVTKJS_CONVERTER=${VTKJS_CONVERTER}
        -DBINARY_PATH=${MeshTest_BINARY_PATH}
        -DVTKJS_OUTPUT_PATH=${PROJECT_SOURCE_DIR}/web/static/vis/${MeshTest_PATH}
        "-DVIS_FILES=${MeshTest_VIS}"
        -DGLOB_MODE=${GLOB_MODE}
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
        --debug-output
        WORKING_DIRECTORY ${MeshTest_SOURCE_PATH}
    )
    set_tests_properties(${TESTER_NAME} PROPERTIES DEPENDS ${TEST_NAME})

endfunction()
