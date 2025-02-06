# cmake-lint: disable=C0103

if(OGS_USE_PIP)
    set(LOCAL_VIRTUALENV_DIR ${PROJECT_BINARY_DIR}/.venv CACHE INTERNAL "")
    set(_venv_bin_dir "bin")
    if(MSVC)
        set(_venv_bin_dir "Scripts")
    endif()
    set(LOCAL_VIRTUALENV_BIN_DIR ${LOCAL_VIRTUALENV_DIR}/${_venv_bin_dir}
        CACHE INTERNAL ""
    )
    set(Python_ROOT_DIR ${LOCAL_VIRTUALENV_DIR})
    set(CMAKE_REQUIRE_FIND_PACKAGE_Python TRUE)
    find_program(UV_TOOL_PATH uv)
    if(UV_TOOL_PATH)
        set(_pip_install_command ${UV_TOOL_PATH} pip install --prefix
                                 ${LOCAL_VIRTUALENV_DIR} CACHE INTERNAL ""
        )
        set(_pip_uninstall_command ${UV_TOOL_PATH} pip uninstall --prefix
                                   ${LOCAL_VIRTUALENV_DIR} CACHE INTERNAL ""
        )
        set(_venv_tool "uv")
    else()
        set(_pip_install_command ${LOCAL_VIRTUALENV_BIN_DIR}/pip install
            CACHE INTERNAL ""
        )
        set(_pip_uninstall_command ${LOCAL_VIRTUALENV_BIN_DIR}/pip uninstall
                                   --yes CACHE INTERNAL ""
        )
        set(_venv_tool "pip")
    endif()
    if(NOT EXISTS ${LOCAL_VIRTUALENV_DIR})
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                -Dpython_version=${ogs.minimum_version.python}...<3.14
                -DUV_TOOL_PATH=${UV_TOOL_PATH} -P
                ${PROJECT_SOURCE_DIR}/scripts/cmake/PythonCreateVirtualEnv.cmake
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR} COMMAND_ECHO STDOUT
                              ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
            RESULT_VARIABLE _return_code
        )
        if(NOT ${_return_code} EQUAL 0)
            message(
                FATAL_ERROR
                    "Creation of Python virtual environment failed!\n"
                    "To disable virtual environments set OGS_USE_PIP=OFF."
            )
        endif()
        unset(_OGS_PYTHON_PACKAGES_SHA1 CACHE)
    endif()
    # Fixes macOS install issues
    execute_process(
        COMMAND ${_pip_install_command} wheel
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    # Create jupytext config
    file(
        WRITE
        ${LOCAL_VIRTUALENV_DIR}/etc/jupyter/labconfig/default_setting_overrides.json
        [=[
{
  "@jupyterlab/docmanager-extension:plugin": {
    "defaultViewers": {
      "markdown": "Jupytext Notebook",
      "myst": "Jupytext Notebook"
    }
  }
}
]=]
    )
else()
    # Prefer unix location over frameworks (Apple-only)
    set(Python_FIND_FRAMEWORK LAST)

    # Prefer more recent Python version
    set(Python_FIND_STRATEGY VERSION)
endif()

set(_python_componets Interpreter Development.Module)
if(NOT OGS_BUILD_WHEEL)
    list(APPEND _python_componets Development.Embed)
endif()

find_package(
    Python ${ogs.minimum_version.python} COMPONENTS ${_python_componets}
    REQUIRED
)

if(OGS_USE_PIP)
    set(Python_SITEARCH_NATIVE ${Python_SITEARCH})
    if(WIN32)
        string(REPLACE "\\" "\\\\" Python_SITEARCH_NATIVE
                       ${Python_SITEARCH_NATIVE}
        )
    endif()
    set(OGS_PYTHON_PACKAGES ""
        CACHE INTERNAL "List of Python packages to be installed via pip."
    )
    set(Python_ROOT_DIR ${LOCAL_VIRTUALENV_DIR})
    if(MSVC)
        set(Python_EXECUTABLE ${Python_ROOT_DIR}/Scripts/python.exe)
    else()
        set(Python_EXECUTABLE ${Python_ROOT_DIR}/bin/python)
    endif()
    if(OGS_BUILD_TESTING)
        # Notebook requirements from Tests/Data
        file(STRINGS Tests/Data/requirements.txt _requirements)
        # \; are not preserved in list operations, substitute via a placeholder
        string(REPLACE "\;" "_semicolon_" _requirements "${_requirements}")
        file(STRINGS Tests/Data/requirements-dev.txt _requirements_dev)
        list(APPEND OGS_PYTHON_PACKAGES ${_requirements} ${_requirements_dev})

        list(APPEND OGS_PYTHON_PACKAGES
             "setuptools" # https://github.com/glenfant/stopit/issues/32
        )
    endif()
endif()

# Sets up a Python virtual environment in the build directory
function(setup_venv)
    # Caches a hash of requested Python packages when they were successfully
    # installed. On subsequent runs compare new hash to cached. If equal do
    # nothing.
    get_property(
        _addtest_python_packages GLOBAL PROPERTY AddTest_PYTHON_PACKAGES
    )
    list(APPEND OGS_PYTHON_PACKAGES ${_addtest_python_packages})
    list(REMOVE_DUPLICATES OGS_PYTHON_PACKAGES)
    list(SORT OGS_PYTHON_PACKAGES)
    string(SHA1 _ogs_python_packages_sha1 "${OGS_PYTHON_PACKAGES}")
    list(LENGTH OGS_PYTHON_PACKAGES OGS_PYTHON_PACKAGES_LENGTH)
    if(NOT ${_ogs_python_packages_sha1} STREQUAL "${_OGS_PYTHON_PACKAGES_SHA1}"
       AND ${OGS_PYTHON_PACKAGES_LENGTH} GREATER 0
    )
        string(REPLACE ";" "\n" REQUIREMENTS_CONTENT "${OGS_PYTHON_PACKAGES}")
        # Revert back _semicolon_ placeholder
        string(REPLACE "_semicolon_" "\\;" REQUIREMENTS_CONTENT
                       "${REQUIREMENTS_CONTENT}"
        )
        file(WRITE ${PROJECT_BINARY_DIR}/requirements.txt
             ${REQUIREMENTS_CONTENT}
        )
        message(
            STATUS
                "Installing Python packages into local virtual environment using ${_venv_tool} ..."
        )
        if(APPLE)
            # CC=/Library/Developer/CommandLineTools/usr/bin/cc and this somehow
            # breaks wheel builds ...
            set(_apple_env ${CMAKE_COMMAND} -E env CC=clang CXX=clang)
        endif()
        execute_process(
            COMMAND ${_apple_env} ${_pip_install_command} -r requirements.txt
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            RESULT_VARIABLE _return_code
            OUTPUT_VARIABLE _out
            ERROR_VARIABLE _err
        )
        if(${_return_code} EQUAL 0)
            set(_OGS_PYTHON_PACKAGES_SHA1 "${_ogs_python_packages_sha1}"
                CACHE INTERNAL ""
            )
            message(STATUS "${_out}")
        else()
            message(
                FATAL_ERROR
                    "Installation of Python packages via pip failed!\n"
                    "To disable pip set OGS_USE_PIP=OFF.\n\n${_out}\n${_err}"
            )
        endif()
        if(DEFINED ENV{CI} AND UNIX AND NOT APPLE)
            set(_pip_gmsh_flags --force-reinstall --pre)
            if(UV_TOOL_PATH)
                set(_pip_gmsh_flags --reinstall --prerelease=allow)
            endif()
            execute_process(
                COMMAND
                    ${_apple_env} ${_pip_install_command} ${_pip_gmsh_flags}
                    --index-url https://gmsh.info/python-packages-dev-nox
                    "gmsh>=4.11"
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                RESULT_VARIABLE _return_code
                OUTPUT_VARIABLE _out
                ERROR_VARIABLE _err
            )
        endif()
        if(${_return_code} EQUAL 0)
            set(_OGS_PYTHON_PACKAGES_SHA1 "${_ogs_python_packages_sha1}"
                CACHE INTERNAL ""
            )
            message(STATUS "${_out}")
        else()
            message(
                FATAL_ERROR
                    "Installation of Python packages via pip failed!\n"
                    "To disable pip set OGS_USE_PIP=OFF.\n\n${_out}\n${_err}"
            )
        endif()
        # Uninstall ogs wheel
        execute_process(
            COMMAND ${_apple_env} ${_pip_uninstall_command} ogs
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    endif()
endfunction()

# Sets up ctest which are dependent on the virtual env, e.g. using ogstools
function(setup_venv_dependent_ctests)
    if(NOT OGS_USE_MPI AND OGS_BUILD_TESTING AND OGS_BUILD_PROCESS_HT)
        execute_process(
            COMMAND
                ${Python_EXECUTABLE}
                ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/generateInvalidMediaForHT.py
            WORKING_DIRECTORY
                ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles
            RESULT_VARIABLE GEN_INVALID_RES
        )
        if(NOT GEN_INVALID_RES EQUAL 0)
            message(SEND_ERROR "generateInvalidMediaForHT.py failed with status ${GEN_INVALID_RES}.")
        endif()
        file(GLOB HT_INVALID_PRJ_FILES
             ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/*.prj
        )
        foreach(ht_invalid_prj_file ${HT_INVALID_PRJ_FILES})
            string(
                REPLACE ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/HT
                        "invalid" ht_invalid_prj_file_short
                        ${ht_invalid_prj_file}
            )
            AddTest(
                NAME HT_${ht_invalid_prj_file_short}
                PATH Parabolic/HT/InvalidProjectFiles
                EXECUTABLE ogs
                EXECUTABLE_ARGS ${ht_invalid_prj_file}
                RUNTIME 1 PROPERTIES WILL_FAIL TRUE
            )
        endforeach()
    endif()
endfunction()
