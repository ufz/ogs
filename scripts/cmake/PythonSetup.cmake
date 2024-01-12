# cmake-lint: disable=C0103

set(_python_version_max "...<3.12")
if(WIN32 AND NOT OGS_BUILD_WHEEL)
    # 3.11 crashes at initialization on Windows.
    set(_python_version_max "...<3.11")
endif()

if(OGS_USE_PIP)
    set(LOCAL_VIRTUALENV_DIR ${PROJECT_BINARY_DIR}/.venv CACHE INTERNAL "")
    set(Python_ROOT_DIR ${LOCAL_VIRTUALENV_DIR})
    set(CMAKE_REQUIRE_FIND_PACKAGE_Python TRUE)
    if(NOT EXISTS ${LOCAL_VIRTUALENV_DIR})
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                -Dpython_version=${ogs.minimum_version.python}${_python_version_max}
                -P
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
    set(_venv_bin_dir "bin")
    if(MSVC)
        set(_venv_bin_dir "Scripts")
    endif()
    set(LOCAL_VIRTUALENV_BIN_DIR ${LOCAL_VIRTUALENV_DIR}/${_venv_bin_dir}
        CACHE INTERNAL ""
    )
    # Fixes macOS install issues
    execute_process(
        COMMAND ${LOCAL_VIRTUALENV_BIN_DIR}/pip install wheel
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
    Python ${ogs.minimum_version.python}${_python_version_max}
    COMPONENTS ${_python_componets} REQUIRED
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
             "snakemake==${ogs.minimum_version.snakemake}"
             "pulp==2.7.0" # https://github.com/snakemake/snakemake/issues/2607
        )
        set(SNAKEMAKE ${LOCAL_VIRTUALENV_BIN_DIR}/snakemake CACHE FILEPATH ""
                                                                  FORCE
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
                "Installing Python packages into local virtual environment..."
        )
        if(APPLE)
            # CC=/Library/Developer/CommandLineTools/usr/bin/cc and this somehow
            # breaks wheel builds ...
            set(_apple_env ${CMAKE_COMMAND} -E env CC=clang CXX=clang)
        endif()
        execute_process(
            COMMAND ${_apple_env} ${LOCAL_VIRTUALENV_BIN_DIR}/pip install -r
                    requirements.txt
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            RESULT_VARIABLE _return_code
            OUTPUT_VARIABLE _out
            ERROR_VARIABLE _err
        )
        if(DEFINED ENV{CI} AND UNIX AND NOT APPLE)
            execute_process(
                COMMAND ${_apple_env} ${LOCAL_VIRTUALENV_BIN_DIR}/pip install
                    --force-reinstall
                    -r ${PROJECT_SOURCE_DIR}/Tests/Data/requirements-gmsh-nox.txt
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
    endif()
endfunction()
