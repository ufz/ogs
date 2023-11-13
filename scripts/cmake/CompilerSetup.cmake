# C++ standard setup
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set compiler helper variables
if((${CMAKE_CXX_COMPILER_ID} MATCHES "Clang") OR (${CMAKE_CXX_COMPILER_ID}
                                                  MATCHES "IntelLLVM")
)
    set(COMPILER_IS_CLANG TRUE CACHE BOOL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(COMPILER_IS_GCC TRUE CACHE BOOL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(COMPILER_IS_INTEL TRUE CACHE BOOL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    set(COMPILER_IS_MSVC TRUE CACHE BOOL "")
endif() # CMAKE_CXX_COMPILER_ID

if(APPLE AND "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm64")
    set(APPLE_ARM TRUE CACHE BOOL "Apple M processors" FORCE)
endif()

# GNU-like compiler
if(COMPILER_IS_GCC OR COMPILER_IS_CLANG OR COMPILER_IS_INTEL)
    # Coloring output
    option(FORCE_COLORED_OUTPUT
           "Always produce ANSI-colored output (GNU/Clang only)." ON
    )
    if(${FORCE_COLORED_OUTPUT})
        if(COMPILER_IS_GCC)
            add_compile_options(-fdiagnostics-color=always)
        elseif(COMPILER_IS_CLANG)
            add_compile_options(-fcolor-diagnostics)
        endif()
    endif()

    # Profiling
    if(OGS_PROFILE)
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
            message(
                STATUS "When using profiling you should set CMAKE_BUILD_TYPE \
                to Release."
            )
        endif()
        set(PROFILE_FLAGS -pg -fno-omit-frame-pointer -O2 -DNDEBUG)
        # clang compiler does not know the following flags
        if(NOT COMPILER_IS_CLANG)
            set(PROFILE_FLAGS
                ${PROFILE_FLAGS} -fno-inline-functions-called-once
                -fno-optimize-sibling-calls
            )
        endif()
        add_compile_options(${PROFILE_ARGS})
    endif()

    if(OGS_CPU_ARCHITECTURE STREQUAL "generic")
        set(CPU_FLAGS -mtune=generic)
    elseif(NOT APPLE_ARM AND OGS_CPU_ARCHITECTURE)
        set(CPU_FLAGS -march=${OGS_CPU_ARCHITECTURE})
    endif()

    if(COMPILER_IS_GCC)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${ogs.minimum_version.gcc})
            message(FATAL_ERROR "Aborting: GCC ${ogs.minimum_version.gcc} is \
                    required! Found version ${CMAKE_CXX_COMPILER_VERSION}."
            )
        endif()
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-fext-numeric-literals>)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 13.1.1
           OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 13.2.1
        )
            # See https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/4597
            add_compile_options(
                $<$<COMPILE_LANGUAGE:CXX>:-Wno-dangling-reference>
                $<$<COMPILE_LANGUAGE:CXX>:-Wno-array-bounds>
                $<$<COMPILE_LANGUAGE:CXX>:-Wno-stringop-overflow>
                $<$<COMPILE_LANGUAGE:CXX>:-Wno-stringop-overread>
            )
        endif()
    endif()

    if(COMPILER_IS_CLANG)
        if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS
               ${ogs.minimum_version.apple_clang}
            )
                message(
                    FATAL_ERROR
                        "Aborting: Apple Clang ${ogs.minimum_version.apple_clang} \
                    is required! Found version ${CMAKE_CXX_COMPILER_VERSION}. Update Xcode!"
                )
            endif()
        else()
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS
               ${ogs.minimum_version.clang}
            )
                message(
                    FATAL_ERROR "Aborting: Clang ${ogs.minimum_version.clang} \
                    is required! Found version ${CMAKE_CXX_COMPILER_VERSION}"
                )
            endif()
        endif()
    endif()

    if(COMPILER_IS_INTEL)
        # Use highest instruction set available on the compilation host
        # processor
        add_compile_options(-xHOST)
    endif()

    # Linker: prefer lld > gold > regular
    foreach(linker lld gold)
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -fuse-ld=${linker} -Wl,--version
            ERROR_QUIET OUTPUT_VARIABLE _linker_version
        )
        if("${_linker_version}" MATCHES "LLD")
            add_link_options(-fuse-ld=lld)
            message(STATUS "Using lld linker. (${_linker_version})")
            break()
        elseif("${_linker_version}" MATCHES "GNU gold")
            add_link_options(-fuse-ld=gold)
            message(STATUS "Using GNU gold linker. (${_linker_version})")
            break()
        endif()
    endforeach()
endif()

if(MSVC)
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS
       ${ogs.minimum_version.msvc.compiler}
    )
        message(
            FATAL_ERROR
                "Aborting: Visual Studio compiler \
            ${ogs.minimum_version.msvc.compiler} is required. Found version \
            ${CMAKE_CXX_COMPILER_VERSION}."
        )
    endif()
    if(${MSVC_TOOLSET_VERSION} LESS ${ogs.minimum_version.msvc.toolset})
        message(
            FATAL_ERROR
                "Aborting: Visual Studio ${ogs.minimum_version.msvc.year} \
            is required! Found Visual Studio with toolset version \
            ${MSVC_TOOLSET_VERSION}. See the following link for version info: \
            https://cmake.org/cmake/help/v3.16/variable/MSVC_TOOLSET_VERSION.html"
        )
    endif()
    if(OGS_CPU_ARCHITECTURE STREQUAL "native")
        set(CPU_FLAGS /favor:blend)
    else()
        set(CPU_FLAGS /favor:${OGS_CPU_ARCHITECTURE})
    endif()
    add_compile_options(
        /wd4290
        /wd4267
        /wd4996
        /bigobj
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_WARNINGS
        -D_CRT_XNONSTDC_NO_WARNINGS
        -D__restrict__=__restrict # this fixes #5
        # This fixes compile errors with std::numeric_limits<T>::min() / max()
        -DNOMINMAX
        -DBOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE # when VC is newer than Boost
        # Disables all warnings coming from include with <>-syntax
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory/
        /experimental:external
        /external:anglebrackets
        /external:W0
    )
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4099")

    # Use Multi-ToolTask scheduler
    if(NOT CMAKE_VS_GLOBALS MATCHES "(^|;)UseMultiToolTask=")
        list(APPEND CMAKE_VS_GLOBALS UseMultiToolTask=true)
    endif()
    if(NOT CMAKE_VS_GLOBALS MATCHES "(^|;)EnforceProcessCountAcrossBuilds=")
        list(APPEND CMAKE_VS_GLOBALS EnforceProcessCountAcrossBuilds=true)
    endif()
endif()

if(PROJECT_IS_TOP_LEVEL)
    include(Sanitizers)
endif()

add_compile_options(
    ${OGS_CXX_FLAGS} # user-given, CMake-option
    ${CPU_FLAGS}
)
