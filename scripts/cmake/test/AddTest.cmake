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
#   WRAPPER <time|mpirun> # optional
#   WRAPPER_ARGS <arguments> # optional
#   TESTER <diff|vtkdiff|vtkdiff-mesh|gmldiff|memcheck|numdiff> # optional
#   TESTER_ARGS <argument> # optional
#   REQUIREMENTS # optional simple boolean expression which has to be true to
#                  enable the test, e.g.
#                  OGS_USE_PETSC AND (FOO OR BAR)
#   PYTHON_PACKAGES package_x=1.2.3 package_y=0.1.x # optional
#   RUNTIME <in seconds> # optional for optimizing ctest duration
#                          values should be taken from envinf job
#   WORKING_DIRECTORY # optional, specify the working directory of the test
#   DISABLED # optional, disables the test
#   PROPERTIES <test properties> # optional
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
        TESTER_ARGS
        REQUIREMENTS
        PYTHON_PACKAGES
        PROPERTIES
        LABELS
    )
    cmake_parse_arguments(
        AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    set(AddTest_SOURCE_PATH "${Data_SOURCE_DIR}/${AddTest_PATH}")
    set(AddTest_BINARY_PATH "${Data_BINARY_DIR}/${AddTest_PATH}")
    set(AddTest_STDOUT_FILE_PATH
        "${AddTest_BINARY_PATH}/${AddTest_NAME}_stdout.txt"
    )

    # set defaults
    if(NOT DEFINED AddTest_EXECUTABLE)
        message(FATAL_ERROR "Test ${AddTest_NAME}: No EXECUTABLE set!")
    endif()
    if(NOT DEFINED AddTest_REQUIREMENTS)
        set(AddTest_REQUIREMENTS TRUE)
    endif()
    if(NOT DEFINED AddTest_RUNTIME)
        set(AddTest_RUNTIME 1)
    elseif(AddTest_RUNTIME GREATER 750)
        # Set a timeout on jobs larger than the default ctest timeout of 1500 (s).
        # The allowed runtime is twice as long as the given RUNTIME parameter.
        math(EXPR timeout "${AddTest_RUNTIME} * 2")
        set(timeout TIMEOUT ${timeout})
    endif()
    if(NOT DEFINED AddTest_WORKING_DIRECTORY)
        set(AddTest_WORKING_DIRECTORY ${AddTest_BINARY_PATH})
    endif()

    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
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

    if(DEFINED OGS_EXCLUDE_CTESTS)
        foreach(regexp ${OGS_EXCLUDE_CTESTS})
            if("${AddTest_NAME}" MATCHES "${regexp}")
                message(
                    STATUS "Disabled by OGS_EXCLUDE_CTESTS: ${AddTest_NAME}"
                )
                return()
            endif()
        endforeach()
    endif()

    # check requirements, disable if not met
    if(NOT TARGET ${AddTest_EXECUTABLE})
        return()
    endif()
    if(${AddTest_REQUIREMENTS})
        message(DEBUG "Enabling test ${AddTest_NAME}.")
    else()
        return()
    endif()

    # --- Implement wrappers ---
    if(AddTest_WRAPPER STREQUAL "time")
        if(TIME_TOOL_PATH)
            set(WRAPPER_COMMAND time)
        else()
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
    if(AddTest_TESTER STREQUAL "vtkdiff-mesh" AND NOT TARGET vtkdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "xdmfdiff" AND NOT TARGET xdmfdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "gmldiff" AND NOT ${Python_Interpreter_FOUND})
        return()
    endif()
    if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "numdiff" AND NOT NUMDIFF_TOOL_PATH)
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

    if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER MATCHES "vtkdiff"
        OR AddTest_TESTER STREQUAL "xdmfdiff") AND NOT AddTest_DIFF_DATA
    )
        message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
        set(TESTER_ARGS "-sbB")
    elseif(AddTest_TESTER MATCHES "vtkdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
    elseif(AddTest_TESTER STREQUAL "xdmfdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:xdmfdiff>)
    elseif(AddTest_TESTER STREQUAL "numdiff")
        set(SELECTED_DIFF_TOOL_PATH ${NUMDIFF_TOOL_PATH})
    endif()

    # -----------
    if(TARGET ${AddTest_EXECUTABLE})
        set(AddTest_EXECUTABLE_PARSED $<TARGET_FILE:${AddTest_EXECUTABLE}>)
    else()
        set(AddTest_EXECUTABLE_PARSED ${AddTest_EXECUTABLE})
    endif()

    # Run the wrapper
    if(DEFINED AddTest_WRAPPER)
        set(AddTest_WRAPPER_STRING "-${AddTest_WRAPPER}")
    endif()
    set(TEST_NAME
        "${AddTest_EXECUTABLE}-${AddTest_NAME}${AddTest_WRAPPER_STRING}"
    )

    # Process placeholders <PATH> <SOURCE_PATH> and <BUILD_PATH>
    string(REPLACE "<PATH>" "${AddTest_PATH}" AddTest_WORKING_DIRECTORY
                   "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<PATH>" "${AddTest_PATH}" AddTest_EXECUTABLE_ARGS
                   "${AddTest_EXECUTABLE_ARGS}"
    )

    string(REPLACE "<SOURCE_PATH>" "${Data_SOURCE_DIR}/${AddTest_PATH}"
                   AddTest_WORKING_DIRECTORY "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<SOURCE_PATH>" "${Data_SOURCE_DIR}/${AddTest_PATH}"
                   AddTest_EXECUTABLE_ARGS "${AddTest_EXECUTABLE_ARGS}"
    )
    string(REPLACE "<BUILD_PATH>" "${Data_BINARY_DIR}/${AddTest_PATH}"
                   AddTest_WORKING_DIRECTORY "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<BUILD_PATH>" "${Data_BINARY_DIR}/${AddTest_PATH}"
                   AddTest_EXECUTABLE_ARGS "${AddTest_EXECUTABLE_ARGS}"
    )

    current_dir_as_list(ProcessLib labels)
    if(DEFINED AddTest_LABELS)
        list(APPEND labels ${AddTest_LABELS})
    else()
        list(APPEND labels default)
    endif()
    if(${AddTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
        list(APPEND labels small)
    else()
        list(APPEND labels large)
    endif()

    if(AddTest_PYTHON_PACKAGES)
        list(APPEND labels python_modules)
        if(OGS_USE_PIP)
            # Info has to be passed by global property because it is not
            # possible to set cache variables from inside a function.
            set_property(
                GLOBAL APPEND PROPERTY AddTest_PYTHON_PACKAGES
                                       ${AddTest_PYTHON_PACKAGES}
            )
        else()
            message(
                STATUS
                    "Warning: Benchmark ${AddTest_NAME} requires these "
                    "Python packages: ${AddTest_PYTHON_PACKAGES}!\n Make sure to "
                    "have them installed in your current Python environment OR "
                    "set OGS_USE_PIP=ON!"
            )
        endif()
    endif()

    _add_test(${TEST_NAME})

    # OpenMP tests for specific processes only. TODO (CL) Once all processes can
    # be assembled OpenMP parallel, the condition should be removed.
    if("${labels}" MATCHES "TH2M|ThermoRichards")
        _add_test(${TEST_NAME}-omp)
        _set_omp_test_properties()
    endif()

    add_dependencies(ctest ${AddTest_EXECUTABLE})
    add_dependencies(ctest-large ${AddTest_EXECUTABLE})

    if(OGS_WRITE_BENCHMARK_COMMANDS)
        string(
            REPLACE
                ";"
                " "
                _cmd
                "${AddTest_WRAPPER} ${AddTest_WRAPPER_ARGS} ${AddTest_EXECUTABLE} ${AddTest_EXECUTABLE_ARGS}"
        )
        string(REPLACE "${Data_BINARY_DIR}/" "" _cmd "${_cmd}")
        string(REPLACE "-o" "" _cmd "${_cmd}")
        # string(STRIP "${_cmd}" _cmd)
        set(_benchmark_run_commands
            "${_benchmark_run_commands} ${AddTest_EXECUTABLE} ; ${_cmd} ; ${TEST_NAME}\n"
            CACHE INTERNAL ""
        )
    endif()

    if(NOT AddTest_TESTER OR OGS_COVERAGE)
        return()
    endif()

    # Run the tester
    _add_test_tester(${TEST_NAME})
    if("${labels}" MATCHES "TH2M|ThermoRichards")
        _add_test_tester(${TEST_NAME}-omp)
    endif()

endfunction()

# Add a ctest and sets properties
macro(_add_test TEST_NAME)

    set(_binary_path ${AddTest_BINARY_PATH})
    if("${TEST_NAME}" MATCHES "-omp")
        set(_binary_path ${_binary_path}-omp)
    endif()

    file(TO_NATIVE_PATH "${_binary_path}" AddTest_BINARY_PATH_NATIVE)
    file(MAKE_DIRECTORY ${_binary_path})
    set(_exe_args ${AddTest_EXECUTABLE_ARGS})
    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(_exe_args -o ${AddTest_BINARY_PATH_NATIVE}
                      ${AddTest_EXECUTABLE_ARGS}
        )
    endif()

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${AddTest_EXECUTABLE_PARSED}
            "-DEXECUTABLE_ARGS=${_exe_args}" # Quoted because
            # passed as list see https://stackoverflow.com/a/33248574/80480
            -DBINARY_PATH=${_binary_path} -DWRAPPER_COMMAND=${WRAPPER_COMMAND}
            "-DWRAPPER_ARGS=${AddTest_WRAPPER_ARGS}"
            -DWORKING_DIRECTORY=${AddTest_WORKING_DIRECTORY}
            -DLOG_FILE=${PROJECT_BINARY_DIR}/logs/${TEST_NAME}.txt -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    if(DEFINED AddTest_DEPENDS)
        if(NOT (TEST ${AddTest_DEPENDS} OR TARGET ${AddTest_DEPENDS}))
            message(
                FATAL_ERROR
                    "AddTest ${TEST_NAME}: dependency ${AddTest_DEPENDS} does not exist!"
            )
        endif()
        set_tests_properties(${TEST_NAME} PROPERTIES DEPENDS ${AddTest_DEPENDS})
    endif()
    if(DEFINED MPI_PROCESSORS)
        set_tests_properties(
            ${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCESSORS}
        )
    endif()

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${AddTest_PROPERTIES}
                   COST
                   ${AddTest_RUNTIME}
                   DISABLED
                   ${AddTest_DISABLED}
                   LABELS
                   "${labels}"
                   ${timeout}
    )
endmacro()

# Sets number of threads, adds label 'omp'
macro(_set_omp_test_properties)
    get_test_property(${TEST_NAME}-omp ENVIRONMENT _environment)
    if(NOT _environment)
        set(_environment "")
    endif()
    set_tests_properties(
        ${TEST_NAME}-omp
        PROPERTIES ENVIRONMENT "OGS_ASM_THREADS=4;${_environment}" PROCESSORS 4
                   LABELS "${labels};omp"
    )
endmacro()

# Adds subsequent tester ctest.
macro(_add_test_tester TEST_NAME)

    unset(TESTER_COMMAND)
    set(_binary_path ${AddTest_BINARY_PATH})
    if("${TEST_NAME}" MATCHES "-omp")
        set(_binary_path ${_binary_path}-omp)
    endif()

    set(TESTER_NAME "${TEST_NAME}-${AddTest_TESTER}")

    if(AddTest_TESTER STREQUAL "diff")
        foreach(FILE ${AddTest_DIFF_DATA})
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(
                APPEND
                TESTER_COMMAND
                "${SELECTED_DIFF_TOOL_PATH} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${_binary_path}/${FILE}"
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
                ${_binary_path}/${VTK_FILE} \
                -a ${NAME_A} -b ${NAME_B} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS}"
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
                    ${_binary_path}/${VTK_FILE} \
                    -a ${NAME_A} -b ${NAME_B} \
                    --abs ${ABS_TOL} --rel ${REL_TOL} \
                    ${TESTER_ARGS} ${AddTest_TESTER_ARGS}"
                    )
                endif()
            endforeach()
        else()
            message(
                FATAL_ERROR
                    "For vtkdiff tester the number of diff data arguments must be a multiple of six."
            )
        endif()
    elseif(AddTest_TESTER STREQUAL "vtkdiff-mesh")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod3 "${DiffDataLength} % 3")
        if(${DiffDataLengthMod3} EQUAL 0)
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOLERANCE)

                list(
                    APPEND
                    TESTER_COMMAND
                    "${SELECTED_DIFF_TOOL_PATH} -m \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${_binary_path}/${VTK_FILE} \
                --abs ${ABS_TOLERANCE}"
                )
            endforeach()
        else()
            message(FATAL_ERROR "The number of diff data arguments must be a
            multiple of three: expected.vtu output.vtu absolute_tolerance."
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
                file(TO_NATIVE_PATH "${Python_EXECUTABLE}" PY_EXE)
                # Dirty hack for Windows Python paths with spaces:
                string(REPLACE "Program Files" "\"Program Files\"" PY_EXE
                               ${PY_EXE}
                )
            else()
                set(PY_EXE ${Python_EXECUTABLE})
            endif()
            list(
                APPEND
                TESTER_COMMAND
                "${PY_EXE} ${PROJECT_SOURCE_DIR}/scripts/test/gmldiff.py \
                --abs ${ABS_TOL} --rel ${REL_TOL} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS}\
                ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${_binary_path}/${GML_FILE}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "memcheck")
        set(TESTER_COMMAND
            "! ${GREP_TOOL_PATH} definitely ${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.txt"
        )
    elseif(AddTest_TESTER STREQUAL "numdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
            list(GET AddTest_DIFF_DATA "${DiffDataIndex}" FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(
                APPEND
                TESTER_COMMAND
                "${SELECTED_DIFF_TOOL_PATH} -a ${ABS_TOL} -r ${REL_TOL}  \
        ${TESTER_ARGS} ${AddTest_TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
        ${_binary_path}/${FILE}"
            )
        endforeach()
    endif()

    add_test(
        NAME ${TESTER_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DSOURCE_PATH=${AddTest_SOURCE_PATH}
            -DSELECTED_DIFF_TOOL_PATH=${SELECTED_DIFF_TOOL_PATH}
            "-DTESTER_COMMAND=${TESTER_COMMAND}" -DBINARY_PATH=${_binary_path}
            -DGLOB_MODE=${GLOB_MODE}
            -DLOG_FILE_BASE=${PROJECT_BINARY_DIR}/logs/${TESTER_NAME} -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
            --debug-output
        WORKING_DIRECTORY ${AddTest_SOURCE_PATH}
    )
    set_tests_properties(
        ${TESTER_NAME} PROPERTIES DEPENDS ${TEST_NAME} DISABLED
                                  ${AddTest_DISABLED} LABELS "tester;${labels}"
    )
endmacro()
