if(NOT DEFINED ENV{GUIX_ENVIRONMENT} AND NOT DEFINED ENV{NIX_BUILD_TOP})
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
set(OGS_BUILD_TESTING OFF CACHE BOOL "" FORCE) # TODO: not yet supported
# TODO: create a newer eigen package:
set(OGS_USE_EIGEN_UNSUPPORTED OFF CACHE BOOL "" FORCE)
set(OGS_INSTALL_DEPENDENCIES OFF CACHE BOOL "" FORCE) # handled by guix
set(OGS_CPU_ARCHITECTURE OFF CACHE BOOL "" FORCE) # enables guix --tune
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
