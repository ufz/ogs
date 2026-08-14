include(FindPackageHandleStandardArgs)

if(TARGET FastFloat::fast_float)
    set(FastFloat_FOUND TRUE)
    return()
endif()

find_package(FastFloat CONFIG QUIET)
if(TARGET FastFloat::fast_float)
    set(FastFloat_FOUND TRUE)
    return()
endif()

find_path(FastFloat_INCLUDE_DIR NAMES fast_float/fast_float.h)

if(FastFloat_INCLUDE_DIR)
    set(_fast_float_version_header
        "${FastFloat_INCLUDE_DIR}/fast_float/float_common.h"
    )
    if(EXISTS "${_fast_float_version_header}")
        file(STRINGS "${_fast_float_version_header}" _fast_float_version_lines
             REGEX "^#define FASTFLOAT_VERSION_(MAJOR|MINOR|PATCH) "
        )
        foreach(_fast_float_version_part MAJOR MINOR PATCH)
            string(
                REGEX MATCH
                      "#define FASTFLOAT_VERSION_${_fast_float_version_part} +([0-9]+)"
                      _fast_float_version_match "${_fast_float_version_lines}"
            )
            set(_fast_float_version_${_fast_float_version_part}
                "${CMAKE_MATCH_1}"
            )
        endforeach()
        set(FastFloat_VERSION
            "${_fast_float_version_MAJOR}.${_fast_float_version_MINOR}.${_fast_float_version_PATCH}"
        )
    endif()
endif()

find_package_handle_standard_args(
    FastFloat
    REQUIRED_VARS FastFloat_INCLUDE_DIR
    VERSION_VAR FastFloat_VERSION
)

if(FastFloat_FOUND AND NOT TARGET FastFloat::fast_float)
    add_library(FastFloat::fast_float INTERFACE IMPORTED)
    target_include_directories(
        FastFloat::fast_float SYSTEM INTERFACE "${FastFloat_INCLUDE_DIR}"
    )
endif()
