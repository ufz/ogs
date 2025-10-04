string(REPLACE ":" ";" PATH_LIST "$ENV{PATH}")

foreach(path_entry IN LISTS PATH_LIST)
    message(STATUS "Checking: ${path_entry}")
    if("${path_entry}" MATCHES "^/gnu/store/[a-z0-9]*-profile/(bin|sbin)$")
        set(MATCHED_GUIX_PATH "${path_entry}")
        message(STATUS "Guix dev environment detected: ${MATCHED_GUIX_PATH}")
        break()
    endif()
endforeach()

if(NOT DEFINED ENV{GUIX_ENVIRONMENT} AND NOT DEFINED ENV{NIX_BUILD_TOP} AND NOT DEFINED MATCHED_GUIX_PATH)
    return()
endif()

if(DEFINED ENV{GUIX_ENVIRONMENT})
    message(STATUS "Guix environment detected: $ENV{GUIX_ENVIRONMENT}")
endif()
if(DEFINED ENV{NIX_BUILD_TOP})
    message(STATUS "Guix build detected: $ENV{NIX_BUILD_TOP}")
    string(REGEX MATCH "^/gnu/store/(.*)-(.*)-(.*)$" _guix_version_match
                 "${CMAKE_INSTALL_PREFIX}"
    )
    if("${OGS_VERSION}" STREQUAL "NO_VERSION")
        set(OGS_VERSION "${CMAKE_MATCH_3}-guix-${CMAKE_MATCH_1}")
        message(STATUS "Using guix-provided OGS_VERSION=${OGS_VERSION}.")
    endif()
endif()

set(GUIX_BUILD ON CACHE BOOL "" FORCE)
if(NOT DEFINED MATCHED_GUIX_PATH)
    set(OGS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
endif()
set(OGS_INSTALL_DEPENDENCIES OFF CACHE BOOL "" FORCE) # handled by guix
set(OGS_CPU_ARCHITECTURE OFF CACHE BOOL "" FORCE) # enables guix --tune
