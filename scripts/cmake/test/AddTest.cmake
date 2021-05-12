#
# AddTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# ~~~
# AddTest(
#   NAME <name of the the test>
#   PATH <working directory> # relative to SourceDir/Tests/Data
#   EXECUTABLE <executable target> # optional, defaults to ogs
#   EXECUTABLE_ARGS <arguments>
#   WRAPPER <time|memcheck|callgrind|mpirun> # optional
#   WRAPPER_ARGS <arguments> # optional
#   TESTER <diff|vtkdiff|gmldiff|memcheck> # optional
#   REQUIREMENTS # optional simple boolean expression which has to be true to
#                  enable the test, e.g.
#                  OGS_USE_PETSC AND (FOO OR BAR)
#   PYTHON_PACKAGES package_x=1.2.3 package_y=0.1.x # optional
#   VIS <vtu output file(s)> # optional for documentation
#   RUNTIME <in seconds> # optional for optimizing ctest duration
#                          values should be taken from envinf job
#   WORKING_DIRECTORY # optional, specify the working directory of the test
#   DISABLED # optional, disables the test
# )
# ~~~
#
# Conditional arguments:
#
# ~~~
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
#
#   gmldiff-tester
#     - DIFF_DATA
#         <gml file> <absolute tolerance> <relative tolerance>
#         Can be given multiple times; the point coordinates in the gml files are
#         compared using the given absolute and relative tolerances.
# ~~~
# cmake-lint: disable=C0103,R0911,R0912,R0915
function(AddTest)

    # parse arguments
    set(options DISABLED)
    set(oneValueArgs
        EXECUTABLE
        PATH
        NAME
        WRAPPER
        TESTER
        ABSTOL
        RELTOL
        RUNTIME
        DEPENDS
        WORKING_DIRECTORY
    )
    set(multiValueArgs
        EXECUTABLE_ARGS
        DATA
        DIFF_DATA
        WRAPPER_ARGS
        REQUIREMENTS
        PYTHON_PACKAGES
        VIS
    )
    cmake_parse_arguments(
        AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    set(AddTest_SOURCE_PATH "${Data_SOURCE_DIR}/${AddTest_PATH}")
    set(AddTest_BINARY_PATH "${Data_BINARY_DIR}/${AddTest_PATH}")
    file(MAKE_DIRECTORY ${AddTest_BINARY_PATH})
    file(TO_NATIVE_PATH "${AddTest_BINARY_PATH}" AddTest_BINARY_PATH_NATIVE)
    set(AddTest_STDOUT_FILE_PATH
        "${AddTest_BINARY_PATH}/${AddTest_NAME}_stdout.log"
    )

    # set defaults
    if(NOT DEFINED AddTest_EXECUTABLE)
        message(FATAL_ERROR "Test ${AddTest_NAME}: No EXECUTABLE set!")
    endif()
    if(NOT DEFINED AddTest_REQUIREMENTS)
        set(AddTest_REQUIREMENTS TRUE)
    endif()
    set(timeout ${ogs.ctest.large_runtime})
    if(DEFINED AddTest_RUNTIME)
        math(EXPR timeout "${AddTest_RUNTIME} * 3")
    else()
        set(AddTest_RUNTIME 1)
    endif()
    if(NOT DEFINED AddTest_WORKING_DIRECTORY)
        set(AddTest_WORKING_DIRECTORY ${AddTest_BINARY_PATH})
    endif()

    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(AddTest_EXECUTABLE_ARGS
            -o ${AddTest_BINARY_PATH_NATIVE}
            ${AddTest_EXECUTABLE_ARGS}
        )
        set(AddTest_WORKING_DIRECTORY ${AddTest_SOURCE_PATH})
    endif()

    if(DEFINED OGS_CTEST_MAX_RUNTIME)
        if(${AddTest_RUNTIME} GREATER ${OGS_CTEST_MAX_RUNTIME})
            return()
        endif()
    endif()
    if(${AddTest_RUNTIME} GREATER ${ogs.ctest.large_runtime})
        string(PREPEND AddTest_NAME "LARGE_")
    endif()

    # --- Implement wrappers ---
    # check if exe is part of build
    if(NOT TARGET ${AddTest_EXECUTABLE})
        set(DISABLED_TESTS_LOG
            "${DISABLED_TESTS_LOG}\nTest exe ${AddTest_EXECUTABLE} not built! Disabling test ${AddTest_NAME}."
            CACHE INTERNAL ""
        )
        return()
    endif()
    # check requirements, disable if not met
    if(${AddTest_REQUIREMENTS})
        message(DEBUG "Enabling test ${AddTest_NAME}.")
    else()
        set(DISABLED_TESTS_LOG
            "${DISABLED_TESTS_LOG}\nRequirement ${AddTest_REQUIREMENTS} not met! Disabling test ${AddTest_NAME}."
            CACHE INTERNAL ""
        )
        return()
    endif()

    if(AddTest_WRAPPER STREQUAL "time")
        if(TIME_TOOL_PATH)
            set(WRAPPER_COMMAND time)
        else()
            set(DISABLED_TESTS_LOG
                "${DISABLED_TESTS_LOG}\nDisabling time wrapper for ${AddTest_NAME} as time exe was not found!"
                CACHE INTERNAL ""
            )
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "memcheck")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND
                "${VALGRIND_TOOL_PATH} --tool=memcheck --log-file=${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log -v --leak-check=full --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xff"
            )
            set(tester memcheck)
        else()
            set(DISABLED_TESTS_LOG
                "${DISABLED_TESTS_LOG}\nDisabling memcheck wrapper for ${AddTest_NAME} as memcheck exe was not found!"
                CACHE INTERNAL ""
            )
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "callgrind")
        if(VALGRIND_TOOL_PATH)
            set(WRAPPER_COMMAND
                "${VALGRIND_TOOL_PATH} --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes"
            )
            unset(tester)
        else()
            set(DISABLED_TESTS_LOG
                "${DISABLED_TESTS_LOG}\nDisabling callgrind wrapper for ${AddTest_NAME} as callgrind exe was not found!"
                CACHE INTERNAL ""
            )
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "mpirun")
        if(MPIRUN_TOOL_PATH)
            if("${HOSTNAME}" MATCHES "frontend.*")
                set(AddTest_WRAPPER_ARGS ${AddTest_WRAPPER_ARGS} --mca
                                         btl_openib_allow_ib 1
                )
            endif()
            set(WRAPPER_COMMAND ${MPIRUN_TOOL_PATH})
            if("${AddTest_WRAPPER_ARGS}" MATCHES "-np;([0-9]*)")
                set(MPI_PROCESSORS ${CMAKE_MATCH_1})
            endif()
        else()
            message(
                STATUS
                    "ERROR: mpirun was not found but is required for ${AddTest_NAME}!"
            )
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
    if(AddTest_TESTER STREQUAL "xdmfdiff" AND NOT TARGET xdmfdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "gmldiff" AND NOT ${Python3_Interpreter_FOUND})
        return()
    endif()
    if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
        return()
    endif()

    if(AddTest_DIFF_DATA)
        string(LENGTH "${AddTest_DIFF_DATA}" DIFF_DATA_LENGTH)
        if(${DIFF_DATA_LENGTH} GREATER 7500)
            message(
                FATAL_ERROR
                    "${AddTest_NAME}: DIFF_DATA to long! Consider using regex-syntax: TODO"
            )
        endif()
    endif()

    if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER STREQUAL "vtkdiff"
        OR AddTest_TESTER STREQUAL "xdmfdiff") AND NOT AddTest_DIFF_DATA
    )
        message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
        set(TESTER_ARGS "-sbB")
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
    elseif(AddTest_TESTER STREQUAL "xdmfdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:xdmfdiff>)
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        foreach(FILE ${AddTest_DIFF_DATA})
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(APPEND TESTER_COMMAND "${SELECTED_DIFF_TOOL_PATH} \
                ${TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${AddTest_BINARY_PATH}/${FILE}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "vtkdiff" OR AddTest_TESTER STREQUAL
                                                "xdmfdiff"
    )
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod4 "${DiffDataLength} % 4")
        math(EXPR DiffDataLengthMod6 "${DiffDataLength} % 6")
        if(${DiffDataLengthMod4} EQUAL 0 AND NOT ${DiffDataLengthMod6} EQUAL 0)
            message(WARNING "DEPRECATED AddTest call with four arguments.\
Use six arguments version of AddTest with absolute and relative tolerances"
            )
            if(NOT AddTest_ABSTOL)
                set(AddTest_ABSTOL 1e-16)
            endif()
            if(NOT AddTest_RELTOL)
                set(AddTest_RELTOL 1e-16)
            endif()
            set(TESTER_ARGS "--abs ${AddTest_ABSTOL} --rel ${AddTest_RELTOL}")
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 4)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)

                list(
                    APPEND
                    TESTER_COMMAND
                    "${SELECTED_DIFF_TOOL_PATH} \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${AddTest_BINARY_PATH}/${VTK_FILE} \
                -a ${NAME_A} -b ${NAME_B} \
                ${TESTER_ARGS}"
                )
            endforeach()
        elseif(${DiffDataLengthMod6} EQUAL 0)
            if(${AddTest_ABSTOL} OR ${AddTest_RELTOL})
                message(
                    FATAL_ERROR
                        "ABSTOL or RELTOL arguments must not be present."
                )
            endif()
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 6)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
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
                    list(APPEND TESTER_COMMAND
                         "${VTK_FILE} ${NAME_A} ${NAME_B} ${ABS_TOL} ${REL_TOL}"
                    )
                    set(GLOB_MODE TRUE)
                else()
                    list(
                        APPEND
                        TESTER_COMMAND
                        "${SELECTED_DIFF_TOOL_PATH} \
                    ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                    ${AddTest_BINARY_PATH}/${VTK_FILE} \
                    -a ${NAME_A} -b ${NAME_B} \
                    --abs ${ABS_TOL} --rel ${REL_TOL} \
                    ${TESTER_ARGS}"
                    )
                endif()
            endforeach()
        else()
            message(
                FATAL_ERROR
                    "For vtkdiff tester the number of diff data arguments must be a multiple of six."
            )
        endif()
    elseif(AddTest_TESTER STREQUAL "gmldiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
            list(GET AddTest_DIFF_DATA "${DiffDataIndex}" GML_FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)

            get_filename_component(FILE_EXPECTED ${GML_FILE} NAME)
            if(WIN32)
                string(REPLACE " " "\\ " PY_EXE ${Python3_EXECUTABLE})
            else()
                set(PY_EXE ${Python3_EXECUTABLE})
            endif()
            list(
                APPEND
                TESTER_COMMAND
                "${PY_EXE} ${PROJECT_SOURCE_DIR}/scripts/test/gmldiff.py \
                --abs ${ABS_TOL} --rel ${REL_TOL} \
                ${TESTER_ARGS} \
                ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${AddTest_BINARY_PATH}/${GML_FILE}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "memcheck")
        set(TESTER_COMMAND
            "! ${GREP_TOOL_PATH} definitely ${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.log"
        )
    endif()

    # -----------
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
    set(TEST_NAME
        "${AddTest_EXECUTABLE}-${AddTest_NAME}${AddTest_WRAPPER_STRING}"
    )
    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${AddTest_EXECUTABLE_PARSED}
            "-DEXECUTABLE_ARGS=${AddTest_EXECUTABLE_ARGS}" # Quoted because
                                                           # passed as list see
                                                           # https://stackoverflow.com/a/33248574/80480
            -DBINARY_PATH=${AddTest_BINARY_PATH}
            -DWRAPPER_COMMAND=${WRAPPER_COMMAND}
            "-DWRAPPER_ARGS=${AddTest_WRAPPER_ARGS}"
            "-DFILES_TO_DELETE=${FILES_TO_DELETE}"
            -DWORKING_DIRECTORY=${AddTest_WORKING_DIRECTORY} -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )
    if(DEFINED AddTest_DEPENDS)
        set_tests_properties(${TEST_NAME} PROPERTIES DEPENDS ${AddTest_DEPENDS})
    endif()
    if(DEFINED MPI_PROCESSORS)
        set_tests_properties(
            ${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCESSORS}
        )
    endif()

    current_dir_as_list(ProcessLib labels)
    if(${AddTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
        list(APPEND labels default)
    else()
        list(APPEND labels large)
    endif()

    set_tests_properties(${TEST_NAME}
        PROPERTIES
            COST ${AddTest_RUNTIME}
            DISABLED ${AddTest_DISABLED}
            LABELS "${labels}"
    )
    # Disabled for the moment, does not work with CI under load
    # if(NOT OGS_COVERAGE)
    #     set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT ${timeout})
    # endif()

    add_dependencies(ctest ${AddTest_EXECUTABLE})
    add_dependencies(ctest-large ${AddTest_EXECUTABLE})

    if(AddTest_PYTHON_PACKAGES)
        if(POETRY)
            file(WRITE ${PROJECT_BINARY_DIR}/tmp_poetry_add.bat
                 "poetry add ${AddTest_PYTHON_PACKAGES}"
            )
            if(WIN32)
                set(EXEC_CMD tmp_poetry_add.bat)
            else()
                set(EXEC_CMD ${BASH_TOOL_PATH} tmp_poetry_add.bat)
            endif()
            execute_process(
                COMMAND ${EXEC_CMD} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            )
        else()
            message(
                STATUS
                    "Warning: Benchmark ${AddTest_NAME} requires these "
                    "Python packages: ${AddTest_PYTHON_PACKAGES}!\n Make sure to "
                    "have them installed in your current Python environment OR "
                    "install the Poetry package manager for Python!"
            )
        endif()
    endif()

    if(NOT AddTest_TESTER OR OGS_COVERAGE)
        return()
    endif()

    # Run the tester
    set(TESTER_NAME "${TEST_NAME}-${AddTest_TESTER}")
    add_test(
        NAME ${TESTER_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DSOURCE_PATH=${AddTest_SOURCE_PATH}
            -DBINARY_PATH=${${AddTest_BINARY_PATH}}
            -DSELECTED_DIFF_TOOL_PATH=${SELECTED_DIFF_TOOL_PATH}
            "-DTESTER_COMMAND=${TESTER_COMMAND}"
            -DVTKJS_CONVERTER=${VTKJS_CONVERTER}
            -DBINARY_PATH=${AddTest_BINARY_PATH}
            -DVTKJS_OUTPUT_PATH=${PROJECT_SOURCE_DIR}/web/static/vis/${AddTest_PATH}
            "-DVIS_FILES=${AddTest_VIS}" -DGLOB_MODE=${GLOB_MODE} -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
            --debug-output
        WORKING_DIRECTORY ${AddTest_SOURCE_PATH}
    )
    set_tests_properties(${TESTER_NAME}
        PROPERTIES
            DEPENDS ${TEST_NAME}
            DISABLED ${AddTest_DISABLED}
            LABELS "tester;${labels}"
    )

endfunction()
