# cmake-lint: disable=C0103

# Sets up ctest that are dependent on the virtual env, e.g. using ogstools
function(setup_venv_dependent_ctests)
    if(NOT OGS_USE_MPI AND OGS_BUILD_TESTING AND OGS_BUILD_PROCESS_HT)
        execute_process(
            COMMAND
                uv run python
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

message(STATUS "┌─ PythonSetup.cmake")
list(APPEND CMAKE_MESSAGE_INDENT "│    ")

set(_python_componets Interpreter Development.Module)
if(NOT OGS_BUILD_WHEEL)
    list(APPEND _python_componets Development.Embed)
endif()

if(OGS_USE_PIP)
    set(LOCAL_VIRTUALENV_DIR ${PROJECT_BINARY_DIR}/.venv CACHE INTERNAL "")
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

        setup_venv_dependent_ctests()
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
      "myst": "Jupytext Notebook"
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
