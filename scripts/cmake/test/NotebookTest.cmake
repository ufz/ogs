#
# NotebookTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# ~~~
# NotebookTest(
#   NAME <name of the the test>
#   NOTEBOOKFILE <path to the notebook file>
#   PYTHON_PACKAGES package_x=1.2.3 package_y=0.1.x # optional, additional
#                                                     Python packages to install
#   RUNTIME <in seconds> # optional for optimizing ctest duration
#                          values should be taken from envinf job
#   LABELS <labelA;labelB;...> # optional, defaults to "default"
#   PROPERTIES <test properties> # optional
#   SKIP_WEB # optional, skips the web site generation
# )
# ~~~
#
# cmake-lint: disable=C0103,R0915,R0912
function(NotebookTest)

    if(NOT OGS_BUILD_CLI OR NOT OGS_BUILD_TESTING OR NOT OGS_USE_PIP)
        return()
    endif()

    set(options DISABLED SKIP_WEB)
    set(oneValueArgs NOTEBOOKFILE RUNTIME)
    set(multiValueArgs PROPERTIES LABELS PYTHON_PACKAGES)
    cmake_parse_arguments(
        NotebookTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    get_filename_component(
        NotebookTest_DIR "${NotebookTest_NOTEBOOKFILE}" DIRECTORY
    )
    get_filename_component(
        NotebookTest_NAME "${NotebookTest_NOTEBOOKFILE}" NAME
    )
    get_filename_component(
        NotebookTest_NAME_WE "${NotebookTest_NOTEBOOKFILE}" NAME_WE
    )

    if(NotebookTest_UNPARSED_ARGUMENTS)
        message(
            FATAL_ERROR
                "Unparsed argument(s) '${NotebookTest_UNPARSED_ARGUMENTS}' to NotebookTest call."
        )
    endif()

    set(NotebookTest_SOURCE_DIR "${Data_SOURCE_DIR}/${NotebookTest_DIR}")
    set(_props "")

    if(NOT DEFINED NotebookTest_RUNTIME)
        set(NotebookTest_RUNTIME 1)
    elseif(NotebookTest_RUNTIME GREATER 750)
        # Set a timeout on jobs larger than the default ctest timeout of 1500
        # (s). The allowed runtime is twice as long as the given RUNTIME
        # parameter.
        math(EXPR timeout "${NotebookTest_RUNTIME} * 2")
        set(timeout TIMEOUT ${timeout})
    endif()

    if(DEFINED OGS_CTEST_MAX_RUNTIME)
        # Skip tests that would require too long to execute.
        if(${NotebookTest_RUNTIME} GREATER ${OGS_CTEST_MAX_RUNTIME})
            return()
        endif()
    endif()

    if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/ProcessLib)
        current_dir_as_list(ProcessLib labels)
    endif()
    list(APPEND labels Notebook python)
    if(DEFINED NotebookTest_LABELS)
        list(APPEND labels ${NotebookTest_LABELS})
    else()
        list(APPEND labels default)
    endif()
    # Notebooks are allowed to run longer than usual benchmarks
    math(EXPR _notebook_large_runtime "10 * ${ogs.ctest.large_runtime}")
    if(${NotebookTest_RUNTIME} LESS_EQUAL ${_notebook_large_runtime})
        list(APPEND labels small)
    else()
        list(APPEND labels large)
        string(APPEND NotebookTest_NAME_WE "-LARGE")
    endif()

    set(NotebookTest_BINARY_DIR "${Data_BINARY_DIR}/${NotebookTest_DIR}")
    file(MAKE_DIRECTORY ${NotebookTest_BINARY_DIR})
    file(TO_NATIVE_PATH "${NotebookTest_BINARY_DIR}"
         NotebookTest_BINARY_DIR_NATIVE
    )

    set(TEST_NAME "nb-${NotebookTest_DIR}/${NotebookTest_NAME_WE}")

    if(OGS_USE_PIP AND DEFINED NotebookTest_PYTHON_PACKAGES)
        list(APPEND labels additional_python_modules)
        foreach(_package ${NotebookTest_PYTHON_PACKAGES})
            list(APPEND _uv_run_args --with ${_package})
        endforeach()
    endif()

    set(_exe_args run ${_uv_run_args} python Notebooks/testrunner.py --out ${Data_BINARY_DIR})
    if(NOT NotebookTest_SKIP_WEB)
        list(APPEND _exe_args --hugo)
        if(DEFINED ENV{CI})
            list(APPEND _exe_args --hugo-out ${PROJECT_BINARY_DIR}/web)
        endif()
    endif()
    list(APPEND _exe_args ${NotebookTest_SOURCE_DIR}/${NotebookTest_NAME})

    isTestCommandExpectedToSucceed(${TEST_NAME} ${NotebookTest_PROPERTIES})
    message(DEBUG "Is test '${TEST_NAME}' expected to succeed? → ${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}")

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} ${CMAKE_COMMAND}
            # TODO: only works if notebook is in a leaf directory
            -DEXECUTABLE=${UV_TOOL_PATH}
            "-DEXECUTABLE_ARGS=${_exe_args}"
            -DWORKING_DIRECTORY=${Data_SOURCE_DIR}
            "-DLOG_ROOT=${PROJECT_BINARY_DIR}/logs"
            "-DLOG_FILE_BASENAME=${NotebookTest_NAME_WE}.txt"
            "-DTEST_COMMAND_IS_EXPECTED_TO_SUCCEED=${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
            -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    list(
        APPEND
        _props
        ENVIRONMENT_MODIFICATION
        PATH=path_list_prepend:$<TARGET_FILE_DIR:ogs>
        ${NotebookTest_PROPERTIES}
    )

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${_props}
                   COST
                   ${NotebookTest_RUNTIME}
                   DISABLED
                   ${NotebookTest_DISABLED}
                   LABELS
                   "${labels}"
                   ${timeout}
                   ENVIRONMENT
                   "CI=1;PYDEVD_DISABLE_FILE_VALIDATION=1;UV_PYTHON=$ENV{UV_PYTHON};UV_PROJECT=$ENV{UV_PROJECT};UV_PROJECT_ENVIRONMENT=$ENV{UV_PROJECT_ENVIRONMENT}"
    )

endfunction()
