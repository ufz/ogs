# C++ standard setup
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# < VS 15.0; macOS: https://github.com/sakra/cotire/issues/139
if(MSVC_VERSION LESS 1910 OR APPLE OR ${CMAKE_CXX_COMPILER} MATCHES "clcache")
    set(OGS_USE_PCH OFF CACHE INTERNAL "")
endif()
if(OGS_USE_PCH)
    include(cotire) # compile time reducer
endif()

if(${CMAKE_CXX_COMPILER} MATCHES "clcache" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "clcache does not cache in Debug config!")
endif()

# Set compiler helper variables
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(COMPILER_IS_CLANG TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(COMPILER_IS_GCC TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(COMPILER_IS_INTEL TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    set(COMPILER_IS_MSVC TRUE CACHE INTERNAL "")
endif() # CMAKE_CXX_COMPILER_ID

if(BUILD_SHARED_LIBS)
    # When static libraries are used in some shared libraries it is required
    # that also the static libraries have position independent code.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

    # Enable Windows DLL support.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

### GNU-like compiler
if(COMPILER_IS_GCC OR COMPILER_IS_CLANG OR COMPILER_IS_INTEL)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Set release compiler flags")
        add_compile_options(-O3)
    elseif(NOT STL_NO_DEBUG)
        # Enable assertions in STL in debug mode.
        add_compile_options(
            -D_GLIBCXX_DEBUG
            -D_GLIBCXX_DEBUG_ASSERT
            -D_GLIBCXX_DEBUG_PEDASSERT
            -D_GLIBCXX_DEBUG_VERIFY
        )
    endif()
    add_compile_options(
        -Wall
        -Wextra
    )

    # Coloring output
    option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." ON)
    if (${FORCE_COLORED_OUTPUT})
        if (COMPILER_IS_GCC)
            add_compile_options (-fdiagnostics-color=always)
        elseif (COMPILER_IS_CLANG)
            add_compile_options (-fcolor-diagnostics)
        endif ()
    endif()

    # Profiling
    if(OGS_PROFILE)
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
            message(STATUS "When using profiling you should set CMAKE_BUILD_TYPE \
                to Release.")
        endif()
        set(PROFILE_FLAGS
            -pg
            -fno-omit-frame-pointer
            -O2
            -DNDEBUG
        )
        # clang compiler does not know the following flags
        if(NOT COMPILER_IS_CLANG)
            set(PROFILE_FLAGS ${PROFILE_FLAGS}
                -fno-inline-functions-called-once
                -fno-optimize-sibling-calls
            )
        endif()
        add_compile_options(${PROFILE_ARGS})
    endif()

    if(OGS_ENABLE_AVX2)
        set(CPU_FLAGS -mavx2 -march=core-avx2)
    elseif(OGS_CPU_ARCHITECTURE STREQUAL "generic")
        set(CPU_FLAGS -mtune=generic)
    else()
        set(CPU_FLAGS -march=${OGS_CPU_ARCHITECTURE})
    endif()

    if(COMPILER_IS_GCC)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.3")
            message(FATAL_ERROR "GCC minimum required version is 7.3! You are \
                using ${CMAKE_CXX_COMPILER_VERSION}.")
        endif()
        add_compile_options(-fext-numeric-literals)
    endif()

    if(COMPILER_IS_CLANG)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.5")
            message(FATAL_ERROR "Aborting: Clang 3.5 is required! Found version \
                ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        include(ClangSanitizer)
    endif()

    if(COMPILER_IS_INTEL)
        # Use highest instruction set available on the compilation host processor
        add_compile_options(-xHOST)
    endif()
endif()

if(MSVC)
    if(OGS_CPU_ARCHITECTURE STREQUAL "native")
        set(CPU_FLAGS /favor:blend)
    else()
        set(CPU_FLAGS /favor:${OGS_CPU_ARCHITECTURE})
    endif()
    if(OGS_ENABLE_AVX2)
        set(CPU_FLAGS ${CPU_FLAGS} /arch:AVX2)
    endif()
    add_compile_options(
        /MP # multi-core compilation
        /W3
        /wd4290 /wd4267 /wd4996
        /bigobj
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_WARNINGS
        -D_CRT_XNONSTDC_NO_WARNINGS
        -D__restrict__=__restrict   # this fixes #5
        # This fixes compile errors with
        # std::numeric_limits<T>::min() / max()
        -DNOMINMAX
        -DBOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE # when VC is newer than Boost
        # Disables all warnings coming from include with <>-syntax
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory/
        /experimental:external /external:anglebrackets /external:W0
    )
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4099")
endif()

add_compile_options(
    ${OGS_CXX_FLAGS} # user-given, CMake-option
    ${CPU_FLAGS}
)
