### Git detection ###
find_package(Git REQUIRED)

if(DEFINED ENV{OGS_VERSION})
    set(OGS_VERSION $ENV{OGS_VERSION})
    message(STATUS "OGS VERSION: ${OGS_VERSION} (set via environment)")
endif()

if(NOT IS_GIT_REPO AND NOT OGS_VERSION)
    execute_process(COMMAND ${GIT_EXECUTABLE} status
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE IS_GIT_REPO
        OUTPUT_QUIET)
    if(IS_GIT_REPO GREATER 0)
        set(IS_GIT_REPO FALSE CACHE INTERNAL "")
        if(DEFINED OGS_VERSION)
            message(WARNING "Using user-provided OGS_VERSION; Submodule setup is skipped!")
        else()
            message(FATAL_ERROR "No git repository found at ${PROJECT_SOURCE_DIR}! "
                "Please use git to obtain the source code! See "
                "https://www.opengeosys.org/docs/devguide/getting-started/get-the-source-code/"
                " OR manually set the OGS_VERSION variable.")
        endif()
    else()
        set(IS_GIT_REPO TRUE CACHE INTERNAL "")
    endif()
endif()

if(IS_GIT_REPO AND NOT OGS_VERSION)
    # Get version info from Git, implementation based on
    # https://github.com/tomtom-international/cpp-dependencies
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --long --dirty --always
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE DESCRIBE_RESULT
        OUTPUT_VARIABLE DESCRIBE_STDOUT
    )
    if(DESCRIBE_RESULT EQUAL 0)
        string(STRIP "${DESCRIBE_STDOUT}" DESCRIBE_STDOUT)
        message(STATUS "Git reported this project's version as '${DESCRIBE_STDOUT}'")
        if(DESCRIBE_STDOUT MATCHES "^(.*)-(dirty)$")
          set(DESCRIBE_DIRTY "${CMAKE_MATCH_2}")
          set(DESCRIBE_STDOUT "${CMAKE_MATCH_1}")
        endif()
        if(DESCRIBE_STDOUT MATCHES "^([0-9a-f]+)$")
          set(DESCRIBE_COMMIT_NAME "${CMAKE_MATCH_1}")
          set(DESCRIBE_STDOUT "")
        elseif(DESCRIBE_STDOUT MATCHES "^(.*)-g([0-9a-f]+)$")
          set(DESCRIBE_COMMIT_NAME "g${CMAKE_MATCH_2}")
          set(DESCRIBE_STDOUT "${CMAKE_MATCH_1}")
        endif()
        if(DESCRIBE_STDOUT MATCHES "^(.*)-([0-9]+)$")
          set(DESCRIBE_COMMIT_COUNT "${CMAKE_MATCH_2}")
          set(DESCRIBE_TAG "${CMAKE_MATCH_1}")
          set(DESCRIBE_STDOUT "")
        endif()

        set(OGS_VERSION ${DESCRIBE_TAG})
        if(DESCRIBE_COMMIT_COUNT GREATER 0)
          set(OGS_VERSION "${OGS_VERSION}-${DESCRIBE_COMMIT_COUNT}-${DESCRIBE_COMMIT_NAME}")
        endif()

        if(DESCRIBE_DIRTY)
          string(TIMESTAMP DESCRIBE_DIRTY_TIMESTAMP "%Y%m%d%H%M%S" UTC)
          set(OGS_VERSION "${OGS_VERSION}.dirty.${DESCRIBE_DIRTY_TIMESTAMP}")
        endif()
        message(STATUS "OGS VERSION: ${OGS_VERSION}")
    else()
        message(WARNING "Git repository contains no tags! Please run: git fetch --tags")
    endif()

    # Get git commit
    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%H
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_SHA1
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_SHA1_SHORT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

### Python setup ###
find_program(POETRY poetry)
if(POETRY)
    configure_file(${PROJECT_SOURCE_DIR}/scripts/python/poetry.in.toml
        ${PROJECT_BINARY_DIR}/poetry.toml COPYONLY)
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/pyproject.toml)
        configure_file(${PROJECT_SOURCE_DIR}/scripts/python/pyproject.in.toml
            ${PROJECT_BINARY_DIR}/pyproject.toml)
    endif()
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/.venv)
        if(MSVC)
            set(CMD_PREFIX cmd /C)
        endif()
        execute_process(
            COMMAND ${CMD_PREFIX} ${POETRY} install
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    endif()
    set(Python3_ROOT_DIR ${PROJECT_BINARY_DIR}/.venv)
    set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/bin/python)
    if(MSVC)
        set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/Scripts/python.exe)
    endif()
endif()

if(OGS_USE_PYTHON)
    find_package(Python3 ${ogs.minimum_version.python} COMPONENTS Interpreter Development REQUIRED)
else()
    find_package(Python3 ${ogs.minimum_version.python} COMPONENTS Interpreter)
endif()
if(POETRY)
    set(Python3_VIRTUALENV_SITEPACKAGES
        ${Python3_ROOT_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)
endif()
