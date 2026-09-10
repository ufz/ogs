# cmake-lint: disable=C0103,C0111,R0912,R0915
include(${PROJECT_SOURCE_DIR}/scripts/cmake/test/TestProperties.cmake)

function(OgsTest)
    if(NOT OGS_BUILD_CLI OR NOT OGS_BUILD_TESTING)
        return()
    endif()

    set(options DISABLED NO_OMP_VARIANT NO_TEST_DEFINITION)
    set(oneValueArgs PROJECTFILE RUNTIME NAME_SUFFIX)
    set(multiValueArgs WRAPPER PROPERTIES LABELS PATCH_FILES EXECUTABLE_ARGS
                       FEATURES
    )
    cmake_parse_arguments(
        OgsTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    get_filename_component(OgsTest_DIR "${OgsTest_PROJECTFILE}" DIRECTORY)
    get_filename_component(OgsTest_NAME "${OgsTest_PROJECTFILE}" NAME)
    get_filename_component(OgsTest_NAME_WE "${OgsTest_PROJECTFILE}" NAME_WE)

    if(OgsTest_UNPARSED_ARGUMENTS)
        message(
            FATAL_ERROR
                "Unparsed argument(s) '${OgsTest_UNPARSED_ARGUMENTS}' to OgsTest call."
        )
    endif()

    ogs_resolve_test_runtime(
        "${OgsTest_RUNTIME}"
        "${ogs.ctest.large_runtime}"
        "${OGS_CTEST_MAX_RUNTIME}"
        OgsTest_RUNTIME
        _timeout
        _is_large
        _skip
    )
    if(_timeout)
        set(timeout TIMEOUT ${_timeout})
    endif()
    if(_skip)
        return()
    endif()
    if(_is_large)
        string(APPEND OgsTest_NAME_WE "-LARGE")
    endif()

    if(DEFINED OGS_EXCLUDE_CTESTS)
        foreach(regexp ${OGS_EXCLUDE_CTESTS})
            if("${OgsTest_NAME}" MATCHES "${regexp}")
                message(
                    STATUS "Disabled by OGS_EXCLUDE_CTESTS: ${OgsTest_NAME}"
                )
                return()
            endif()
        endforeach()
    endif()

    set(OgsTest_SOURCE_DIR "${Data_SOURCE_DIR}/${OgsTest_DIR}")
    set(TEST_NAME "ogs-${OgsTest_DIR}/${OgsTest_NAME_WE}")
    set(_processors 1)
    # Add wrapper postfix (-mpi for mpirun).
    if(OgsTest_WRAPPER)
        list(GET OgsTest_WRAPPER 0 WRAPPER)
        if(WRAPPER STREQUAL "mpirun")
            set(TEST_NAME "${TEST_NAME}-mpi")
            list(APPEND OgsTest_WRAPPER --bind-to none)
            if("${OgsTest_WRAPPER}" MATCHES "-np?;([0-9]*)")
                set(_processors ${CMAKE_MATCH_1})
            else()
                message(
                    FATAL_ERROR
                        "${TEST_NAME}: mpirun-wrapper requires -np argument"
                )
            endif()
        endif()
    endif()

    if(DEFINED OgsTest_NAME_SUFFIX)
        set(TEST_NAME "${TEST_NAME}-${OgsTest_NAME_SUFFIX}")
    elseif(OgsTest_PATCH_FILES)
        # Append short hash of the patch files to the test name.
        string(SHA1 _patches_hash "${OgsTest_PATCH_FILES}")
        string(SUBSTRING "${_patches_hash}" 0 4 _short_patches_hash)
        set(TEST_NAME "${TEST_NAME}_${_short_patches_hash}")
        message(DEBUG
                "Test name is already defined. New test name: ${TEST_NAME}"
        )
    endif()

    if(OgsTest_NO_TEST_DEFINITION)
        set(_exe_args ${OgsTest_SOURCE_DIR}/${OgsTest_NAME})
    else()
        set(_exe_args -r ${OgsTest_SOURCE_DIR}
                      ${OgsTest_SOURCE_DIR}/${OgsTest_NAME}
        )
    endif()
    foreach(_patch ${OgsTest_PATCH_FILES})
        list(APPEND _exe_args -p ${OgsTest_SOURCE_DIR}/${_patch})
    endforeach()
    if(OgsTest_EXECUTABLE_ARGS)
        list(APPEND _exe_args ${OgsTest_EXECUTABLE_ARGS})
    endif()

    current_dir_as_list(ProcessLib labels)
    if(OgsTest_LABELS)
        list(APPEND labels ${OgsTest_LABELS})
    else()
        list(APPEND labels default)
    endif()

    if(${OgsTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
        list(APPEND labels small)
    else()
        list(APPEND labels large)
    endif()

    list(APPEND labels ${OgsTest_FEATURES})
    foreach(feature IN LISTS OgsTest_FEATURES)
        if(feature STREQUAL "petsc-mumps")
            # https://gitlab.opengeosys.org/ogs/ogs/-/commit/ff2e3b1024a777a230efb3646890e30ee74004c6
            set(OgsTest_NO_OMP_VARIANT TRUE)
            if(NOT (OGS_USE_PETSC AND OGS_PETSC_HAVE_MUMPS))
                set(OgsTest_DISABLED TRUE)
            endif()
        else()
            message(FATAL_ERROR "Unknown OgsTest feature '${feature}'.")
        endif()
    endforeach()

    set(_has_omp_variant FALSE)
    list(JOIN OGS_OPENMP_PARALLEL_ASM_PROCESSES ";|;"
         match_parallel_asm_processes
    )
    # OpenMP tests for specific processes only. TODO (CL) Once all processes can
    # be assembled OpenMP parallel, the condition should be removed.
    if(NOT OgsTest_NO_OMP_VARIANT AND ";${labels};" MATCHES
                                      ";${match_parallel_asm_processes};"
    )
        set(_has_omp_variant TRUE)
    endif()

    set(_add_non_omp_variant TRUE)
    if(NOT OGS_ENABLE_NON_OMP_TEST_VARIANTS AND _has_omp_variant)
        set(_add_non_omp_variant FALSE)
    endif()

    if(_add_non_omp_variant)
        _ogs_add_test(${TEST_NAME})
    endif()

    if(_has_omp_variant)
        _ogs_add_test(${TEST_NAME}-omp)
        ogs_set_omp_test_properties(
            "${TEST_NAME}-omp" "${_processors}" "${labels}"
        )
    endif()
endfunction()

# Adds a ctest and sets properties
macro(_ogs_add_test TEST_NAME)
    # TEST_NAME is unique, shortened hash added to the working directory of the
    # test to prevent race conditions.
    set(_unique_string "${TEST_NAME}")
    string(SHA1 _unique_hash "${_unique_string}")
    string(SUBSTRING "${_unique_hash}" 0 8 _short_hash)

    # Create unique directory name
    set(OgsTest_BINARY_DIR "${Data_BINARY_DIR}/${OgsTest_DIR}_${_short_hash}")
    file(MAKE_DIRECTORY ${OgsTest_BINARY_DIR})
    file(TO_NATIVE_PATH "${OgsTest_BINARY_DIR}" OgsTest_BINARY_DIR_NATIVE)
    string(REPLACE "/" "_" TEST_NAME_UNDERSCORE ${TEST_NAME})

    isTestCommandExpectedToSucceed(${TEST_NAME} ${OgsTest_PROPERTIES})
    message(
        DEBUG
        "Is test '${TEST_NAME}' expected to succeed? → ${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
    )

    set(_ogs_exe $<TARGET_FILE:ogs>)
    if(OGS_BUILD_WHEEL)
        # When testing the installed wheel assume executable is in PATH from
        # venv.
        set(_ogs_exe ogs)
    endif()
    set(_diff_tool_environment VTKDIFF_EXE=$<TARGET_FILE:vtkdiff>)
    if(TARGET xdmfdiff)
        list(APPEND _diff_tool_environment XDMFDIFF_EXE=$<TARGET_FILE:xdmfdiff>)
    endif()

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${_ogs_exe}
            "-DEXECUTABLE_ARGS=${_exe_args}"
            "-DWRAPPER_COMMAND=${OgsTest_WRAPPER}"
            -DWORKING_DIRECTORY=${OgsTest_BINARY_DIR}
            "-DLOG_FILE_BASENAME=${TEST_NAME_UNDERSCORE}.txt"
            "-DLOG_ROOT=${PROJECT_BINARY_DIR}/logs"
            "-DTEST_COMMAND_IS_EXPECTED_TO_SUCCEED=${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
            -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${OgsTest_PROPERTIES}
                   ENVIRONMENT
                   "${_diff_tool_environment}"
                   COST
                   ${OgsTest_RUNTIME}
                   DISABLED
                   ${OgsTest_DISABLED}
                   LABELS
                   "${labels}"
                   PROCESSORS
                   ${_processors}
                   WORKING_DIRECTORY
                   ${OgsTest_BINARY_DIR}
                   ${timeout}
    )
endmacro()
