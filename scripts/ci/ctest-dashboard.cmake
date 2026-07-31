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
        "\"${CMAKE_COMMAND}\" -E chdir \"${CTEST_SOURCE_DIRECTORY}\" \"${CMAKE_COMMAND}\" --preset=$ENV{CMAKE_PRESET} --log-level=VERBOSE -Wno-dev"
    )
    if(NOT "$ENV{CMAKE_ARGS}" STREQUAL "")
        string(APPEND CTEST_CONFIGURE_COMMAND " $ENV{CMAKE_ARGS}")
    endif()

    ctest_start(Experimental GROUP "${_ctest_group}")
    ctest_configure(
        BUILD "${CTEST_BINARY_DIRECTORY}" SOURCE "${CTEST_SOURCE_DIRECTORY}"
        RETURN_VALUE _configure_result
    )
    ctest_submit(PARTS Configure RETURN_VALUE _submit_result)
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

if(NOT "$ENV{CTEST_DASHBOARD_PHASE}" STREQUAL "build")
    message(
        FATAL_ERROR
            "CTEST_DASHBOARD_PHASE must be set to 'configure' or 'build'."
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
ctest_submit(PARTS Build RETURN_VALUE _submit_result)
if(_build_result)
    message(FATAL_ERROR "CTest build failed with exit code ${_build_result}.")
endif()
if(_submit_result)
    message(
        FATAL_ERROR
            "Submitting CTest build results failed with exit code ${_submit_result}."
    )
endif()
