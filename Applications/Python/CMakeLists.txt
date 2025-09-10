# wheel: Install into Python module root dir (enables 'import ogs.OGSSimulator')
set(_py_install_location ogs)
set(_py_build_location ogs)
if(NOT OGS_BUILD_WHEEL)
    set(_py_install_location
        "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/ogs"
    )
    set(_py_build_location "${PROJECT_BINARY_DIR}/site-packages/ogs")
endif()

if(SKBUILD AND NOT OGS_BUILD_WHEEL)
    message(
        FATAL_ERROR
            "Invalid configuration! Don't use --config-settings=cmake.args!"
            "Use --config-settings=cmake.define instead or supply '--preset wheel` to cmake.args"
    )
endif()

add_subdirectory(ogs)
add_subdirectory(ogs.OGSSimulator)
add_subdirectory(ogs.OGSMesh)
add_subdirectory(ogs.mpl)
if(OGS_BUILD_WHEEL)
    add_subdirectory(ogs.callbacks)
endif()

set_target_properties(
    OGSSimulator OGSMesh mpl PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_py_build_location}
)
file(
    COPY ogs/.
    DESTINATION ${_py_build_location}
    # for guix development container with included (read-only) source code:
    NO_SOURCE_PERMISSIONS # defaults to 644
    PATTERN "__pycache__" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
)

if(OGS_USE_PIP)
    # Fails with libGitInfoLib.so: undefined symbol: __asan_report_load8
    if(NOT ENABLE_ASAN)
        if (OGS_USE_MPI)
            set(OGS_PYTEST_ARGS "-m;not ogs_needs_serial_build")
        endif()
        add_test(
            NAME pytest
            COMMAND
                ${CMAKE_COMMAND} -DEXECUTABLE=${UV_TOOL_PATH}
                "-DEXECUTABLE_ARGS=run;pytest;-ra;${OGS_PYTEST_ARGS};-c;${PROJECT_SOURCE_DIR}/pyproject.toml" # Quoted
                                                                            # because
                # passed as list see https://stackoverflow.com/a/33248574/80480
                -DBINARY_PATH=${_binary_path}
                -DWORKING_DIRECTORY=${PROJECT_SOURCE_DIR}
                "-DLOG_ROOT=${PROJECT_BINARY_DIR}/logs"
                "-DLOG_FILE_BASENAME=pytest.txt"
                "-DTEST_COMMAND_IS_EXPECTED_TO_SUCCEED=TRUE"
                "-DOGS_CI_ALWAYS_SAVE_LOG_FILE_TO_ARTIFACTS=TRUE"
                -P
                ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
        )
        set_tests_properties(pytest PROPERTIES
            LABELS "default;python"
            COST 10
            ENVIRONMENT
                "CI=1;PYDEVD_DISABLE_FILE_VALIDATION=1;UV_PYTHON=$ENV{UV_PYTHON};UV_PROJECT=$ENV{UV_PROJECT};UV_PROJECT_ENVIRONMENT=$ENV{UV_PROJECT_ENVIRONMENT}")
    endif()
endif()
