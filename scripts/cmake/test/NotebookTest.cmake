# cmake-lint: disable=C0103,R0915
function(NotebookTest)

    if(NOT OGS_BUILD_CLI OR NOT OGS_BUILD_TESTING OR NOT OGS_USE_PIP)
        return()
    endif()
    set(options DISABLED)
    set(oneValueArgs NOTEBOOKFILE RUNTIME)
    set(multiValueArgs WRAPPER RESOURCE_LOCK)
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

    set(timeout ${ogs.ctest.large_runtime})
    if(DEFINED NotebookTest_RUNTIME)
        math(EXPR timeout "${NotebookTest_RUNTIME} * 3")
    else()
        set(NotebookTest_RUNTIME 1)
    endif()

    if(DEFINED OGS_CTEST_MAX_RUNTIME)
        # Skip tests that would require too long to execute.
        if(${NotebookTest_RUNTIME} GREATER ${OGS_CTEST_MAX_RUNTIME})
            return()
        endif()
    endif()
    if(${NotebookTest_RUNTIME} GREATER ${ogs.ctest.large_runtime})
        string(APPEND NotebookTest_NAME_WE "-LARGE")
    endif()

    set(NotebookTest_SOURCE_DIR "${Data_SOURCE_DIR}/${NotebookTest_DIR}")
    set(NotebookTest_BINARY_DIR "${Data_BINARY_DIR}/${NotebookTest_DIR}")
    file(MAKE_DIRECTORY ${NotebookTest_BINARY_DIR})
    file(TO_NATIVE_PATH "${NotebookTest_BINARY_DIR}"
         NotebookTest_BINARY_DIR_NATIVE
    )

    set(TEST_NAME "nb-${NotebookTest_DIR}/${NotebookTest_NAME_WE}")

    set(_exe_args Notebooks/testrunner.py --hugo --out ${Data_BINARY_DIR})
    if(DEFINED ENV{CI})
        list(APPEND _exe_args --hugo-out ${PROJECT_BINARY_DIR}/web)
    endif()
    list(APPEND _exe_args ${NotebookTest_SOURCE_DIR}/${NotebookTest_NAME})

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND}
            -E env PYVISTA_HEADLESS=1
            ${CMAKE_COMMAND}
            # TODO: only works if notebook is in a leaf directory
            # -DFILES_TO_DELETE=${Data_BINARY_DIR}/${NotebookTest_DIR}
            -DEXECUTABLE=${Python_EXECUTABLE}
            "-DEXECUTABLE_ARGS=${_exe_args}"
            -DWORKING_DIRECTORY=${Data_SOURCE_DIR} -DCAT_LOG=TRUE -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/OgsTestWrapper.cmake
    )

    current_dir_as_list(ProcessLib labels)
    list(APPEND labels Notebook)
    if(${NotebookTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
        list(APPEND labels default)
    else()
        list(APPEND labels large)
    endif()

    set(_prop_env ENVIRONMENT_MODIFICATION
                  PATH=path_list_prepend:$<TARGET_FILE_DIR:ogs>
    )
    if(DEFINED NotebookTest_RESOURCE_LOCK)
        set_tests_properties(
            ${TEST_NAME}
            PROPERTIES RESOURCE_LOCK ${NotebookTest_RESOURCE_LOCK}
        )
    endif()

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${_prop_env}
                   COST
                   ${NotebookTest_RUNTIME}
                   DISABLED
                   ${NotebookTest_DISABLED}
                   LABELS
                   "${labels}"
    )

endfunction()
