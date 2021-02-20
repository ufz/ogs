### Git detection ###
find_package(Git REQUIRED)

if(DEFINED ENV{OGS_VERSION})
    set(OGS_VERSION $ENV{OGS_VERSION})
    message(STATUS "OGS VERSION: ${OGS_VERSION} (set via environment)")
elseif(DEFINED OGS_VERSION)
    message(STATUS "Using user-provided OGS_VERSION=${OGS_VERSION}")
endif()

if(NOT IS_GIT_REPO)
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

if(IS_GIT_REPO)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE OGS_GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
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
          set(OGS_VERSION "${OGS_VERSION}.dirty")
          if(DEFINED ENV{CI})
            string(TIMESTAMP DESCRIBE_DIRTY_TIMESTAMP "%Y%m%d%H%M%S" UTC)
            set(OGS_VERSION "${OGS_VERSION}.dirty.${DESCRIBE_DIRTY_TIMESTAMP}")
          endif()
        endif()
        message(STATUS "OGS VERSION: ${OGS_VERSION} (reported by git)")
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

if(MSVC)
    set(CMD_COMMAND "cmd" "/c" CACHE INTERNAL "")
endif()
