# cmake-lint: disable=C0103
if(OGS_USE_POETRY)
    find_program(POETRY poetry)
    if(POETRY)
        configure_file(${PROJECT_SOURCE_DIR}/scripts/python/poetry.in.toml
            ${PROJECT_BINARY_DIR}/poetry.toml COPYONLY)
        if(NOT EXISTS ${PROJECT_BINARY_DIR}/pyproject.toml)
            configure_file(${PROJECT_SOURCE_DIR}/scripts/python/pyproject.in.toml
                ${PROJECT_BINARY_DIR}/pyproject.toml)
        endif()
        if(NOT EXISTS ${PROJECT_BINARY_DIR}/.venv)
            execute_process(COMMAND ${_CMD_COMMAND} poetry install
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
        endif()
        set(Python3_ROOT_DIR ${PROJECT_BINARY_DIR}/.venv)
        set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/bin/python)
        if(MSVC)
            set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/Scripts/python.exe)
        endif()
    endif()
endif()

if(OGS_USE_PYTHON)
    find_package(Python3 ${ogs.minimum_version.python} COMPONENTS Interpreter Development REQUIRED)
else()
    find_package(Python3 ${ogs.minimum_version.python} COMPONENTS Interpreter)
endif()
if(POETRY)
    if(MSVC)
        file(TO_NATIVE_PATH "${Python3_ROOT_DIR}/Lib/site-packages"
            Python3_VIRTUALENV_SITEPACKAGES)
        string(REPLACE "\\" "\\\\" Python3_VIRTUALENV_SITEPACKAGES
            ${Python3_VIRTUALENV_SITEPACKAGES})
    else()
        set(Python3_VIRTUALENV_SITEPACKAGES
            ${Python3_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)
    endif()
endif()

set(LOCAL_VIRTUALENV_BIN_DIRS
    ${PROJECT_BINARY_DIR}/.venv/bin
    ${PROJECT_BINARY_DIR}/.venv/Scripts
    CACHE INTERNAL ""
)

if(POETRY)
    if(OGS_BUILD_TESTING)
        list(APPEND PYTHON_PACKAGES snakemake=${ogs.minimum_version.snakemake})
    endif()
    execute_process(COMMAND ${_CMD_COMMAND} poetry add ${PYTHON_PACKAGES}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif()
