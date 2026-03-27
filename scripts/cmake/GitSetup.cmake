# cmake-lint: disable=W0106
set(_ogs_version_user_provided FALSE)
if(DEFINED ENV{OGS_VERSION})
    set(OGS_VERSION $ENV{OGS_VERSION})
    set(_ogs_version_user_provided TRUE)
    message(
        STATUS
            "Using user-provided OGS_VERSION=${OGS_VERSION} (set via environment)."
    )
elseif(DEFINED OGS_VERSION)
    set(_ogs_version_user_provided TRUE)
    message(STATUS "Using user-provided OGS_VERSION=${OGS_VERSION}.")
endif()

# Git detection
find_package(Git QUIET)
if(Git_FOUND)
    # Git for Windows can append a packaging suffix that CMake's version
    # comparison should ignore.
    string(REPLACE ".windows.1" "" GIT_VERSION_STRING ${GIT_VERSION_STRING})
elseif(NOT _ogs_version_user_provided)
    message(
        FATAL_ERROR
            "Git was not found. Either install Git or provide OGS_VERSION as an environment or CMake variable."
    )
else()
    message(
        STATUS
            "Git was not found. Using user-provided OGS_VERSION without Git metadata."
    )
endif()

if(Git_FOUND AND NOT DEFINED _is_git_repo)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _GIT_STATUS
        OUTPUT_VARIABLE _GIT_IS_INSIDE_WORK_TREE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(_GIT_STATUS EQUAL 0 AND _GIT_IS_INSIDE_WORK_TREE STREQUAL "true")
        set(_is_git_repo TRUE)
    else()
        set(_is_git_repo FALSE)
        if(NOT _ogs_version_user_provided)
            if(DEFINED ENV{CI})
                message(
                    FATAL_ERROR
                        "No git repository found at ${PROJECT_SOURCE_DIR}! "
                        "Please use git to obtain the source code OR manually set the OGS_VERSION variable."
                )
            else()
                set(OGS_VERSION "NO_VERSION")
                message(
                    WARNING
                        "No git repository found at ${PROJECT_SOURCE_DIR}. OGS_VERSION is set to NO_VERSION."
                )
            endif()
        endif()
    endif()
endif()

if(Git_FOUND AND _is_git_repo)
    if(DEFINED ENV{CI_COMMIT_BRANCH})
        set(OGS_GIT_BRANCH $ENV{CI_COMMIT_BRANCH})
    else()
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            RESULT_VARIABLE _GIT_BRANCH_RESULT
            OUTPUT_VARIABLE OGS_GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(NOT _GIT_BRANCH_RESULT EQUAL 0)
            unset(OGS_GIT_BRANCH)
            message(WARNING "Could not determine git branch name.")
        endif()
    endif()

    # Get git commit
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _GIT_SHA1_RESULT
        OUTPUT_VARIABLE GIT_SHA1
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT _GIT_SHA1_RESULT EQUAL 0)
        unset(GIT_SHA1)
        message(WARNING "Could not determine git commit SHA1.")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short=8 HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _GIT_SHA1_SHORT_RESULT
        OUTPUT_VARIABLE GIT_SHA1_SHORT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT _GIT_SHA1_SHORT_RESULT EQUAL 0)
        unset(GIT_SHA1_SHORT)
        message(WARNING "Could not determine short git commit SHA1.")
    endif()
endif()

if(Git_FOUND AND _is_git_repo AND NOT _ogs_version_user_provided)
    # Get version info from Git, implementation based on
    # https://github.com/tomtom-international/cpp-dependencies
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --long --dirty --always
                --abbrev=8
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _describe_result
        OUTPUT_VARIABLE _describe_stdout
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(_describe_result EQUAL 0)
        unset(DESCRIBE_DIRTY)
        unset(DESCRIBE_COMMIT_NAME)
        unset(DESCRIBE_COMMIT_COUNT)
        unset(DESCRIBE_TAG)

        if(_describe_stdout MATCHES "^(.*)-(dirty)$")
            set(DESCRIBE_DIRTY "${CMAKE_MATCH_2}")
            set(_describe_stdout "${CMAKE_MATCH_1}")
        endif()
        if(_describe_stdout MATCHES "^([0-9a-f]+)$")
            set(DESCRIBE_COMMIT_NAME "g${CMAKE_MATCH_1}")
            set(_describe_stdout "")
        elseif(_describe_stdout MATCHES "^(.*)-g([0-9a-f]+)$")
            set(DESCRIBE_COMMIT_NAME "g${CMAKE_MATCH_2}")
            set(_describe_stdout "${CMAKE_MATCH_1}")
        endif()
        if(_describe_stdout MATCHES "^(.*)-([0-9]+)$")
            set(DESCRIBE_COMMIT_COUNT "${CMAKE_MATCH_2}")
            set(DESCRIBE_TAG "${CMAKE_MATCH_1}")
            set(_describe_stdout "")
        endif()

        if(DEFINED DESCRIBE_TAG)
            set(OGS_VERSION "${DESCRIBE_TAG}")
            if(DESCRIBE_COMMIT_COUNT GREATER 0)
                set(OGS_VERSION
                    "${OGS_VERSION}-${DESCRIBE_COMMIT_COUNT}-${DESCRIBE_COMMIT_NAME}"
                )
            endif()
        else()
            set(OGS_VERSION "${DESCRIBE_COMMIT_NAME}")
            message(
                WARNING
                    "Git repository contains no tags. Using commit hash as OGS_VERSION."
            )
        endif()

        if(DESCRIBE_DIRTY)
            if(DEFINED ENV{CI})
                string(TIMESTAMP DESCRIBE_DIRTY_TIMESTAMP "%Y%m%d%H%M%S" UTC)
                set(OGS_VERSION "${OGS_VERSION}.dirty.${DESCRIBE_DIRTY_TIMESTAMP}")
            else()
                set(OGS_VERSION "${OGS_VERSION}.dirty")
            endif()
        endif()
        message(STATUS "OGS VERSION: ${OGS_VERSION} (reported by git)")
    else()
        if(DEFINED GIT_SHA1_SHORT)
            set(OGS_VERSION "g${GIT_SHA1_SHORT}")
        else()
            set(OGS_VERSION "NO_VERSION")
        endif()
        message(
            WARNING
                "Could not determine OGS_VERSION from git. Using fallback OGS_VERSION=${OGS_VERSION}. Please run: git fetch --tags"
        )
    endif()
endif()
