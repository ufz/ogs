# cmake-lint: disable=C0103
include(${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTest.cmake)

message(STATUS "┌─ PythonSetup.cmake")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")

set(_python_componets Interpreter Development.Module)
# manylinux_x image used in cibuildwheel on Linux does not contain
# the python library.
if(NOT (LINUX AND DEFINED ENV{CIBUILDWHEEL}))
    list(APPEND _python_componets Development.Embed)
endif()

if(OGS_USE_PIP)
    set(LOCAL_VIRTUALENV_DIR ${PROJECT_BINARY_DIR}/.venv CACHE INTERNAL "")
    # The test environment from Tests/Data is used, not the top-level environment
    # which is used for building wheels only.
    set(ENV{UV_PROJECT} ${PROJECT_SOURCE_DIR}/Tests/Data)
    set(ENV{UV_PROJECT_ENVIRONMENT} ${LOCAL_VIRTUALENV_DIR})
    set(ENV{UV_FROZEN} 1)

    set(_venv_bin_dir "bin")
    if(MSVC)
        set(_venv_bin_dir "Scripts")
    endif()
    set(LOCAL_VIRTUALENV_BIN_DIR ${LOCAL_VIRTUALENV_DIR}/${_venv_bin_dir}
        CACHE INTERNAL ""
    )
    set(Python_ROOT_DIR ${LOCAL_VIRTUALENV_DIR})
    find_program(UV_TOOL_PATH uv REQUIRED)

    # Prefer more recent Python version
    set(Python_FIND_STRATEGY VERSION)
    # Prefer unix location over frameworks (Apple-only)
    set(Python_FIND_FRAMEWORK LAST)

    if(NOT EXISTS ${LOCAL_VIRTUALENV_DIR})
        # Don't use venv
        set(Python_FIND_VIRTUALENV STANDARD)

        find_package(Python ${ogs.minimum_version.python}...<3.14
            COMPONENTS ${_python_componets} REQUIRED)

        set(ENV{UV_PYTHON} ${Python_EXECUTABLE})
        execute_process(
            COMMAND ${UV_TOOL_PATH} venv
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
        execute_process(
            COMMAND ${UV_TOOL_PATH} sync
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            RESULT_VARIABLE _return_code
        )
        if(NOT ${_return_code} EQUAL 0)
            message(
                FATAL_ERROR
                    "Installation of Python packages via uv failed!\n"
                    "To disable uv set OGS_USE_PIP=OFF."
            )
        endif()
    else()
        find_package(Python ${ogs.minimum_version.python}...<3.14
            COMPONENTS ${_python_componets} REQUIRED)
        set(ENV{UV_PYTHON} ${Python_EXECUTABLE})
    endif()

    # Create jupytext config
    file(
        WRITE
        ${LOCAL_VIRTUALENV_DIR}/etc/jupyter/labconfig/default_setting_overrides.json
        [=[
{
  "@jupyterlab/docmanager-extension:plugin": {
    "defaultViewers": {
      "markdown": "Jupytext Notebook",
      "myst": "Jupytext Notebook",
      "python": "Jupytext Notebook"
    }
  }
}
]=]
    )
else()
    find_package(Python ${ogs.minimum_version.python}...<3.14
        COMPONENTS ${_python_componets} REQUIRED)
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "└─ End PythonSetup.cmake")
