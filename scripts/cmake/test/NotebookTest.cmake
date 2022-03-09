# cmake-lint: disable=C0103,R0915
function(NotebookTest)

    if(NOT OGS_BUILD_CLI OR NOT OGS_BUILD_TESTING OR NOT OGS_USE_PIP)
        return()
    endif()
    set(options DISABLED)
    set(oneValueArgs NOTEBOOKFILE RUNTIME)
    set(multiValueArgs WRAPPER)
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

    set(_exe_args Notebooks/testrunner.py --out ${Data_BINARY_DIR}
                  ${NotebookTest_SOURCE_DIR}/${NotebookTest_NAME}
    )

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${Python3_EXECUTABLE}
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

    if(MSVC AND ${CMAKE_VERSION} VERSION_LESS 3.22)
        # ENVIRONMENT_MODIFICATION parameter of set_tests_properties() is
        # required to correctly set the PATH environment variable on Windows.
        message(
            WARNING "Notebook tests are disabled on Windows when CMake < 3.22!"
        )
        return()
    endif()

    if(${CMAKE_VERSION} VERSION_LESS 3.22)
        # This branch applies to *nix only.
        set(_prop_env ENVIRONMENT PATH=$<TARGET_FILE_DIR:ogs>:$ENV{PATH})
    else()
        set(_prop_env ENVIRONMENT_MODIFICATION
                      PATH=path_list_prepend:$<TARGET_FILE_DIR:ogs>
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
