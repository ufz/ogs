# Implementation from https://stackoverflow.com/questions/10205986
#
# Captures not-yet cached CMake variables. On first CMake run via cmake-cli this
# works as expected. Once the variables are cached this will not work anymore
# (and is therefore skipped).
#
# When running CMake -D.. passed args can be retrieved with by the
# CACHE_VARIABLES CMake property. On sub-sequent CMake runs it is no longer
# possible to differentiate between variables already cached by option()-calls
# and variables passed to the CMake call with -D..
#
# A (cached) map data structure would solve this. Tried the following map
# implementations without luck: - https://github.com/toeb/cmakepp (not cached) -
# https://github.com/j3lamp/mcl (did not work at all)

if(EXISTS ${PROJECT_BINARY_DIR}/CMakeCache.txt)
    return()
endif()

get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(cache_var ${CACHE_VARS})
    get_property(CACHE_VAR_HELPSTRING CACHE ${cache_var} PROPERTY HELPSTRING)
    if(CACHE_VAR_HELPSTRING STREQUAL
       "No help, variable specified on the command line."
    )
        get_property(CACHE_VAR_TYPE CACHE ${cache_var} PROPERTY TYPE)
        if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
            set(CACHE_VAR_TYPE)
        else()
            set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
        endif()
        file(TO_CMAKE_PATH "${${cache_var}}" cache_value)
        set(CMAKE_ARGS
            "${CMAKE_ARGS} -D${cache_var}${CACHE_VAR_TYPE}=\"${cache_value}\""
        )
    endif()
endforeach()

set(CMAKE_ARGS "${CMAKE_ARGS} -DCMAKE_BUILD_TYPE=\"${CMAKE_BUILD_TYPE}\""
    CACHE STRING ""
)
string(REPLACE "\"" "\\\"" CMAKE_ARGS_ESCAPED ${CMAKE_ARGS})
set(CMAKE_ARGS_ESCAPED "${CMAKE_ARGS_ESCAPED}" CACHE STRING "")
file(WRITE ${PROJECT_BINARY_DIR}/cmake-args "${CMAKE_ARGS}\n")
