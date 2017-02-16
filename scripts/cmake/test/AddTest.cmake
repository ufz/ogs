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
#   WRAPPER <time|memcheck|callgrind> # optional
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
#       # the given file is compared to a file with the same name from Tests/Data
#
#   numdiff-tester
#     - DIFF_DATA <list of files to numdiff>
#       # the given file is compared to a file with the same name from Tests/Data
#
#   vtkdiff-tester
#     - DIFF_DATA <vtk file> <data array a name> <data array b name>
#       # the given data arrays in the vtk file are compared
#

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
    if (NOT AddTest_ABSTOL)
        set (AddTest_ABSTOL 1e-16)
    endif()
    if (NOT AddTest_RELTOL)
        set (AddTest_RELTOL 1e-16)
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
    if(AddTest_WRAPPER STREQUAL "time" AND NOT TIME_TOOL_PATH)
        return()
    endif()
    if(AddTest_WRAPPER STREQUAL "memcheck" AND NOT VALGRIND_TOOL_PATH)
        return()
    endif()
    if(AddTest_WRAPPER STREQUAL "callgrind" AND NOT VALGRIND_TOOL_PATH)
        return()
    endif()
    if(AddTest_WRAPPER STREQUAL "mpirun" AND NOT MPIRUN_TOOL_PATH)
        return()
    endif()

    if(AddTest_WRAPPER STREQUAL "time")
        set(WRAPPER_COMMAND time)
    elseif(AddTest_WRAPPER STREQUAL "memcheck" AND VALGRIND_TOOL_PATH)
        set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff")
        set(tester memcheck)
    elseif(AddTest_WRAPPER STREQUAL "callgrind" AND VALGRIND_TOOL_PATH)
        set(WRAPPER_COMMAND "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes")
        unset(tester)
    elseif(AddTest_WRAPPER STREQUAL "mpirun")
        set(WRAPPER_COMMAND ${MPIRUN_TOOL_PATH})
    endif()

    # --- Implement testers ---
    # check requirements, disable if not met
    if(AddTest_TESTER STREQUAL "diff" AND NOT DIFF_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "numdiff" AND NOT NUMDIFF_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "vtkdiff" AND NOT TARGET vtkdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
        return()
    endif()

    if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "numdiff" OR AddTest_TESTER STREQUAL "vtkdiff") AND NOT AddTest_DIFF_DATA)
        message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
        set(TESTER_ARGS "-sbB")
    elseif(AddTest_TESTER STREQUAL "numdiff")
        set(SELECTED_DIFF_TOOL_PATH ${NUMDIFF_TOOL_PATH})
        set(TESTER_ARGS "--statistics --absolute-tolerance=${AddTest_ABSTOL} --relative-tolerance=${AddTest_RELTOL}")
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
        set(TESTER_ARGS "--abs ${AddTest_ABSTOL} --rel ${AddTest_RELTOL}")
    endif()

    if(AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "numdiff")
        foreach(FILE ${AddTest_DIFF_DATA})
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} \
                ${TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${AddTest_BINARY_PATH}/${FILE}")
        endforeach()
        string(REPLACE ";" " && " TESTER_COMMAND "${TESTER_COMMAND}")
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod4 "${DiffDataLength} % 4")
        if (NOT ${DiffDataLengthMod4} EQUAL 0)
            message(FATAL_ERROR "For vtkdiff tester the number of diff data arguments must be a multiple of four.")
        endif()

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

        string(REPLACE ";" " && " TESTER_COMMAND "${TESTER_COMMAND}")
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
        "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}"
        -Dcase_path=${AddTest_SOURCE_PATH}
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
        -DTESTER_COMMAND=${TESTER_COMMAND}
        -DVTKJS_CONVERTER=${VTKJS_CONVERTER}
        -DBINARY_PATH=${AddTest_BINARY_PATH}
        -DVTKJS_OUTPUT_PATH=${CMAKE_SOURCE_DIR}/web/static/vis/${AddTest_PATH}
        "-DVIS_FILES=${AddTest_VIS}"
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
    )
    set_tests_properties(${TESTER_NAME} PROPERTIES DEPENDS ${TEST_NAME})

endfunction()
