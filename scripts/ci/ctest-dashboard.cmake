set(CTEST_SOURCE_DIRECTORY "$ENV{CI_PROJECT_DIR}")
get_filename_component(
    CTEST_SOURCE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}" ABSOLUTE
)
get_filename_component(
    CTEST_BINARY_DIRECTORY "$ENV{CTEST_BINARY_DIRECTORY}" ABSOLUTE BASE_DIR
    "${CTEST_SOURCE_DIRECTORY}"
)

include("${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake")
set(CTEST_USE_LAUNCHERS ON)
if(DEFINED BUILDNAME)
    set(CTEST_BUILD_NAME "${BUILDNAME}")
endif()
if(NOT "$ENV{CTEST_SITE}" STREQUAL "")
    set(CTEST_SITE "$ENV{CTEST_SITE}")
else()
    cmake_host_system_information(RESULT CTEST_SITE QUERY HOSTNAME)
endif()

# The CI templates load the generated build environment between these phases.
if("$ENV{CTEST_DASHBOARD_PHASE}" STREQUAL "configure")
    if(DEFINED ENV{CTEST_GROUP} AND NOT "$ENV{CTEST_GROUP}" STREQUAL "")
        set(_ctest_group "$ENV{CTEST_GROUP}")
    else()
        set(_ctest_group Experimental)
    endif()

    set(CTEST_CONFIGURE_COMMAND
        "\"${CMAKE_COMMAND}\" -E chdir \"${CTEST_SOURCE_DIRECTORY}\" \"${CMAKE_COMMAND}\" --preset=$ENV{CMAKE_PRESET} -B \"${CTEST_BINARY_DIRECTORY}\" --log-level=VERBOSE -Wno-dev"
    )
    if(NOT "$ENV{CMAKE_ARGS}" STREQUAL "")
        string(APPEND CTEST_CONFIGURE_COMMAND " $ENV{CMAKE_ARGS}")
    endif()

    ctest_start(Experimental GROUP "${_ctest_group}")
    ctest_configure(
        BUILD "${CTEST_BINARY_DIRECTORY}" SOURCE "${CTEST_SOURCE_DIRECTORY}"
        RETURN_VALUE _configure_result
    )
    if(NOT "$ENV{CTEST_SUBMIT}" STREQUAL "false")
        ctest_submit(PARTS Configure RETURN_VALUE _submit_result)
    else()
        set(_submit_result 0)
    endif()
    if(_configure_result)
        message(
            FATAL_ERROR
                "CTest configure failed with exit code ${_configure_result}."
        )
    endif()
    if(_submit_result)
        message(
            FATAL_ERROR
                "Submitting CTest configure results failed with exit code ${_submit_result}."
        )
    endif()
    return()
endif()

if("$ENV{CTEST_DASHBOARD_PHASE}" STREQUAL "test")
    if(DEFINED ENV{CTEST_GROUP} AND NOT "$ENV{CTEST_GROUP}" STREQUAL "")
        set(_ctest_group "$ENV{CTEST_GROUP}")
    else()
        set(_ctest_group Experimental)
    endif()

    if(EXISTS "${CTEST_BINARY_DIRECTORY}/Testing/TAG")
        ctest_start(APPEND)
    else()
        ctest_start(Experimental GROUP "${_ctest_group}")
    endif()
    set(_ctest_test_command
        ${CMAKE_CTEST_COMMAND} --test-dir "${CTEST_BINARY_DIRECTORY}" -M
        Experimental --group "${_ctest_group}" --no-tests=error -T Test
        --output-junit Tests/ctest.xml
    )
    if(NOT "$ENV{CTEST_ARGS}" STREQUAL "")
        separate_arguments(_ctest_args NATIVE_COMMAND "$ENV{CTEST_ARGS}")
        list(APPEND _ctest_test_command ${_ctest_args})
    elseif(NOT "$ENV{CTEST_PRESET}" STREQUAL "")
        list(APPEND _ctest_test_command --preset "$ENV{CTEST_PRESET}")
    endif()
    if(NOT "$ENV{CTEST_EXTRA_ARGS}" STREQUAL "")
        separate_arguments(
            _ctest_extra_args NATIVE_COMMAND "$ENV{CTEST_EXTRA_ARGS}"
        )
        list(APPEND _ctest_test_command ${_ctest_extra_args})
    endif()
    if(NOT "$ENV{BUILD_CTEST_LARGE}" STREQUAL "true")
        list(APPEND _ctest_test_command -LE large)
    endif()
    if("$ENV{CI_MERGE_REQUEST_LABELS}" MATCHES ".*web[ ]only.*")
        list(APPEND _ctest_test_command -R nb-)
    endif()
    set(_test_timeout)
    if(NOT "$ENV{CTEST_TIMEOUT}" STREQUAL "")
        math(EXPR _test_timeout "$ENV{CTEST_TIMEOUT} * 60")
    endif()

    if(_test_timeout)
        execute_process(
            COMMAND ${_ctest_test_command}
            WORKING_DIRECTORY "${CTEST_SOURCE_DIRECTORY}"
            TIMEOUT ${_test_timeout}
            RESULT_VARIABLE _test_result
        )
    else()
        execute_process(
            COMMAND ${_ctest_test_command}
            WORKING_DIRECTORY "${CTEST_SOURCE_DIRECTORY}"
            RESULT_VARIABLE _test_result
        )
    endif()
    if(_test_result)
        message(FATAL_ERROR "CTest test failed with exit code ${_test_result}.")
    endif()
    if(NOT "$ENV{CTEST_SUBMIT}" STREQUAL "false")
        set(_submit_parts Test)
        set(_buildinfo_file
            "${CTEST_BINARY_DIRECTORY}/Testing/Notes/buildinfo.txt"
        )
        if(EXISTS "${_buildinfo_file}")
            set(CTEST_NOTES_FILES "${_buildinfo_file}")
            list(APPEND _submit_parts Notes)
        endif()
        ctest_submit(PARTS ${_submit_parts} RETURN_VALUE _submit_result)
        if(_submit_result)
            message(
                FATAL_ERROR
                    "Submitting CTest test results failed with exit code ${_submit_result}."
            )
        endif()
    endif()
    return()
endif()

if(NOT "$ENV{CTEST_DASHBOARD_PHASE}" STREQUAL "build")
    message(
        FATAL_ERROR
            "CTEST_DASHBOARD_PHASE must be set to 'configure', 'build', or 'test'."
    )
endif()

ctest_start(APPEND)

set(_build_command_prefix)
if(NOT "$ENV{BUILD_CMD_PREFIX}" STREQUAL "")
    separate_arguments(
        _build_command_prefix NATIVE_COMMAND "$ENV{BUILD_CMD_PREFIX}"
    )
endif()
set(_build_command_without_target
    ${_build_command_prefix} "${CMAKE_COMMAND}" -E chdir
    "${CTEST_SOURCE_DIRECTORY}" "${CMAKE_COMMAND}" --build
    "--preset=$ENV{CMAKE_PRESET}"
)

set(CTEST_BUILD_COMMAND "$ENV{BUILD_CMD_PREFIX}")
if(NOT CTEST_BUILD_COMMAND STREQUAL "")
    string(APPEND CTEST_BUILD_COMMAND " ")
endif()
string(
    APPEND
    CTEST_BUILD_COMMAND
    "\"${CMAKE_COMMAND}\" -E chdir \"${CTEST_SOURCE_DIRECTORY}\" \"${CMAKE_COMMAND}\" --build --preset=$ENV{CMAKE_PRESET}"
)
if("$ENV{BUILD_PACKAGE}" STREQUAL "true")
    string(APPEND CTEST_BUILD_COMMAND " --target package")
endif()
if("$ENV{CMAKE_PRESET}" MATCHES "msvc")
    string(APPEND CTEST_BUILD_COMMAND " -- /m")
endif()

if(DEFINED ENV{ADDITIONAL_TARGETS_PRE} AND NOT "$ENV{ADDITIONAL_TARGETS_PRE}"
                                           STREQUAL ""
)
    separate_arguments(
        _additional_targets NATIVE_COMMAND "$ENV{ADDITIONAL_TARGETS_PRE}"
    )
    foreach(additional_target IN LISTS _additional_targets)
        execute_process(
            COMMAND ${_build_command_without_target} --target
                    "${additional_target}" RESULT_VARIABLE _target_result
        )
        if(_target_result)
            message(
                FATAL_ERROR
                    "Additional target ${additional_target} failed with exit code ${_target_result}."
            )
        endif()
    endforeach()
endif()

ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE _build_result)
if(NOT "$ENV{CTEST_SUBMIT}" STREQUAL "false")
    ctest_submit(PARTS Build RETURN_VALUE _submit_result)
else()
    set(_submit_result 0)
endif()
if(_build_result)
    message(FATAL_ERROR "CTest build failed with exit code ${_build_result}.")
endif()
if(_submit_result)
    message(
        FATAL_ERROR
            "Submitting CTest build results failed with exit code ${_submit_result}."
    )
endif()
