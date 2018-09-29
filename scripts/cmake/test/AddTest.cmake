#
# AddTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# AddTest(
#   NAME <name of the the test>
#   PATH <working directory> # relative to SourceDir/Tests/Data
#   EXECUTABLE <executable target> # optional, defaults to ogs
#   EXECUTABLE_ARGS <arguments>
#   WRAPPER <time|memcheck|callgrind|mpirun> # optional
#   WRAPPER_ARGS <arguments> # optional
#   TESTER <diff|vtkdiff|memcheck> # optional
#   REQUIREMENTS # optional simple boolean expression which has to be true to
#                  enable the test, e.g.
#                  OGS_USE_PETSC AND (OGS_USE_EIGEN OR OGS_USE_LIS)
#   VIS <vtu output file(s)> # optional for documentation
# )
#
# Conditional arguments:
#
#   diff-tester
#     - DIFF_DATA <list of files to diff>
#         the given file is compared to a file with the same name from Tests/Data
#
#   vtkdiff-tester
#     - DIFF_DATA
#         <vtk file a> <vtk file b> <data array a name> <data array b name> <absolute tolerance> <relative tolerance>
#         Can be given multiple times; the given data arrays in the vtk files are
#         compared using the given absolute and relative tolerances.
#       OR
#     - DIFF_DATA
#         GLOB <globbing expression, e.g. xyz*.vtu> <data array a name> <data array b name> <absolute tolerance> <relative tolerance>
#         Searches for all matching files in the working directory (PATH).
#         Matched files are then compared against files with the same name in
#         the benchmark output directory.

function (AddTest)
    if(NOT OGS_BUILD_TESTS)
        return()
    endif()
    # parse arguments
    set(options NONE)
    set(oneValueArgs EXECUTABLE PATH NAME WRAPPER TESTER ABSTOL RELTOL)
    set(multiValueArgs EXECUTABLE_ARGS DATA DIFF_DATA WRAPPER_ARGS REQUIREMENTS VIS)
    cmake_parse_arguments(AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})


    set(AddTest_SOURCE_PATH "${Data_SOURCE_DIR}/${AddTest_PATH}")
    set(AddTest_BINARY_PATH "${Data_BINARY_DIR}/${AddTest_PATH}")
    file(MAKE_DIRECTORY ${AddTest_BINARY_PATH})
    file(TO_NATIVE_PATH "${AddTest_BINARY_PATH}" AddTest_BINARY_PATH_NATIVE)
    set(AddTest_STDOUT_FILE_PATH "${AddTest_BINARY_PATH}/${AddTest_NAME}_stdout.log")

    # set defaults
    if(NOT AddTest_EXECUTABLE)
        set(AddTest_EXECUTABLE ogs)
    endif()
    if (NOT AddTest_REQUIREMENTS)
        set (AddTest_REQUIREMENTS TRUE)
    endif()

    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(AddTest_EXECUTABLE_ARGS -o ${AddTest_BINARY_PATH_NATIVE} ${AddTest_EXECUTABLE_ARGS})
    endif()

    # --- Implement wrappers ---
    # check requirements, disable if not met
    if(${AddTest_REQUIREMENTS})
        # message(STATUS "Enabling test ${AddTest_NAME}.")
    else()
        message(STATUS "Requirement ${AddTest_REQUIREMENTS} not met! Disabling test ${AddTest_NAME}.")
        return()
    endif()

    if(AddTest_WRAPPER STREQUAL "time")
        if(TIME_TOOL_PATH)
            set(WRAPPER_COMMAND time)
        else()
            message(STATUS "WARNING: Disabling time wrapper for ${AddTest_NAME} as time exe was not found!")
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "memcheck")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
            set(tester memcheck)
        else()
            message(STATUS "WARNING: Disabling memcheck wrapper for ${AddTest_NAME} as memcheck exe was not found!")
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "callgrind")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
            unset(tester)
        else()
            message(STATUS "WARNING: Disabling callgrind wrapper for ${AddTest_NAME} as callgrind exe was not found!")
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "mpirun")
        if(MPIRUN_TOOL_PATH)
            set(WRAPPER_COMMAND ${MPIRUN_TOOL_PATH})
        else()
            message(STATUS "ERROR: mpirun was not found but is required for ${AddTest_NAME}!")
            return()
        endif()
    endif()

    # --- Implement testers ---
    # check requirements, disable if not met
    if(AddTest_TESTER STREQUAL "diff" AND NOT DIFF_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "vtkdiff" AND NOT TARGET vtkdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
        return()
    endif()

    if(AddTest_DIFF_DATA)
        string(LENGTH "${AddTest_DIFF_DATA}" DIFF_DATA_LENGTH)
        if(${DIFF_DATA_LENGTH} GREATER 7500)
            message(FATAL_ERROR "${AddTest_NAME}: DIFF_DATA to long! Consider using regex-syntax: TODO")
        endif()
    endif()

    if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "vtkdiff") AND NOT AddTest_DIFF_DATA)
        message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
        set(TESTER_ARGS "-sbB")
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        foreach(FILE ${AddTest_DIFF_DATA})
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} \
                ${TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${AddTest_BINARY_PATH}/${FILE}")
        endforeach()
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod4 "${DiffDataLength} % 4")
        math(EXPR DiffDataLengthMod6 "${DiffDataLength} % 6")
        if (${DiffDataLengthMod4} EQUAL 0 AND NOT ${DiffDataLengthMod6} EQUAL 0)
            message(WARNING "DEPRECATED AddTest call with four arguments.\
Use six arguments version of AddTest with absolute and relative tolerances")
            if (NOT AddTest_ABSTOL)
                set (AddTest_ABSTOL 1e-16)
            endif()
            if (NOT AddTest_RELTOL)
                set (AddTest_RELTOL 1e-16)
            endif()
            set(TESTER_ARGS "--abs ${AddTest_ABSTOL} --rel ${AddTest_RELTOL}")
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 4)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}" REFERENCE_VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)

                list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${AddTest_BINARY_PATH}/${VTK_FILE} \
                -a ${NAME_A} -b ${NAME_B} \
                ${TESTER_ARGS}")
            endforeach()
        elseif (${DiffDataLengthMod6} EQUAL 0)
            if (${AddTest_ABSTOL} OR ${AddTest_RELTOL})
                message(FATAL_ERROR "ABSTOL or RELTOL arguments must not be present.")
            endif()
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 6)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}" REFERENCE_VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+4")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+5")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)

                if("${REFERENCE_VTK_FILE}" STREQUAL "GLOB")
                    list(APPEND TESTER_COMMAND "${VTK_FILE} ${NAME_A} ${NAME_B} ${ABS_TOL} ${REL_TOL}")
                    set(GLOB_MODE TRUE)
                else()
                    list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} \
                    ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                    ${AddTest_BINARY_PATH}/${VTK_FILE} \
                    -a ${NAME_A} -b ${NAME_B} \
                    --abs ${ABS_TOL} --rel ${REL_TOL} \
                    ${TESTER_ARGS}")
                endif()
            endforeach()
        else ()
            message(FATAL_ERROR "For vtkdiff tester the number of diff data arguments must be a multiple of six.")
        endif()
    elseif(tester STREQUAL "memcheck")
        set(TESTER_COMMAND "! ${GREP_TOOL_PATH} definitely ${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log")
    endif()

    ## -----------
    if(TARGET ${AddTest_EXECUTABLE})
        set(AddTest_EXECUTABLE_PARSED $<TARGET_FILE:${AddTest_EXECUTABLE}>)
    else()
        set(AddTest_EXECUTABLE_PARSED ${AddTest_EXECUTABLE})
    endif()

    set(FILES_TO_DELETE "")
    list(APPEND FILES_TO_DELETE "${AddTest_STDOUT_FILE_PATH}")
    foreach(ITEM ${AddTest_DIFF_DATA})
        if(ITEM MATCHES "^.*\.(vtu|vtk)$")
            list(APPEND FILES_TO_DELETE "${ITEM}")
        endif()
    endforeach()

    # Run the wrapper
    if(DEFINED AddTest_WRAPPER)
        set(AddTest_WRAPPER_STRING "-${AddTest_WRAPPER}")
    endif()
    set(TEST_NAME "${AddTest_EXECUTABLE}-${AddTest_NAME}${AddTest_WRAPPER_STRING}")
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
        -DEXECUTABLE=${AddTest_EXECUTABLE_PARSED}
        "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}" # Quoted because passed as list
        -Dcase_path=${AddTest_SOURCE_PATH}             # see https://stackoverflow.com/a/33248574/80480
        -DBINARY_PATH=${AddTest_BINARY_PATH}
        -DWRAPPER_COMMAND=${WRAPPER_COMMAND}
        "-DWRAPPER_ARGS=${AddTest_WRAPPER_ARGS}"
        "-DFILES_TO_DELETE=${FILES_TO_DELETE}"
        -DSTDOUT_FILE_PATH=${AddTest_STDOUT_FILE_PATH}
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    if(TARGET ${AddTest_EXECUTABLE})
        add_dependencies(ctest ${AddTest_EXECUTABLE})
        add_dependencies(ctest-large ${AddTest_EXECUTABLE})
    endif()

    if(NOT AddTest_TESTER OR OGS_COVERAGE)
        return()
    endif()

    # Run the tester
    set(TESTER_NAME "${TEST_NAME}-${AddTest_TESTER}")
    add_test(
        NAME ${TESTER_NAME}
        COMMAND ${CMAKE_COMMAND}
        -Dcase_path=${AddTest_SOURCE_PATH}
        -DBINARY_PATH=${${AddTest_BINARY_PATH}}
        -DSELECTED_DIFF_TOOL_PATH=${SELECTED_DIFF_TOOL_PATH}
        "-DTESTER_COMMAND=${TESTER_COMMAND}"
        -DVTKJS_CONVERTER=${VTKJS_CONVERTER}
        -DBINARY_PATH=${AddTest_BINARY_PATH}
        -DVTKJS_OUTPUT_PATH=${PROJECT_SOURCE_DIR}/web/static/vis/${AddTest_PATH}
        "-DVIS_FILES=${AddTest_VIS}"
        -DGLOB_MODE=${GLOB_MODE}
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
        --debug-output
        WORKING_DIRECTORY ${AddTest_SOURCE_PATH}
    )
    set_tests_properties(${TESTER_NAME} PROPERTIES DEPENDS ${TEST_NAME})

endfunction()
