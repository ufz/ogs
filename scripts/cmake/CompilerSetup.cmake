include(ResetConfigurations)        # To Debug, Release, RelWithDebInfo
include(SetDefaultBuildType)
include(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Debug)
include(MSVCMultipleProcessCompile) # /MP switch (multi processor) for VS
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# C++ standard setup
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC_VERSION LESS 1910) # < VS 15.0
    set(OGS_USE_PCH OFF CACHE INTERNAL "")
endif()
if(OGS_USE_PCH)
    include(cotire) # compile time reducer
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

if(OGS_ENABLE_AVX2)
    set(CPU_FLAGS "-mavx2 -march=core-avx2")
elseif(OGS_CPU_ARCHITECTURE STREQUAL "generic")
    set(CPU_FLAGS "-mtune=generic")
else()
    set(CPU_FLAGS "-march=${OGS_CPU_ARCHITECTURE}")
endif()

### GNU C/CXX compiler
if(COMPILER_IS_GCC)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
        message(FATAL_ERROR "GCC minimum required version is 4.9! You are \
            using ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Set GCC release flags")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
    else()
        # Enable assertions in STL in debug mode.
        if(NOT STL_NO_DEBUG)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
                -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_ASSERT \
                -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY")
        endif()
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPU_FLAGS} -Wno-deprecated -Wall \
        -Wextra -fext-numeric-literals")
endif() # COMPILER_IS_GCC

### Clang
if(COMPILER_IS_CLANG)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.5")
        message(FATAL_ERROR "Aborting: Clang 3.5 is required! Found version \
            ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPU_FLAGS} -Wall \
        -Wno-c++98-compat-pedantic")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        # Enable assertions in STL in debug mode.
        if(NOT STL_NO_DEBUG)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
                -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_ASSERT \
                -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY")
        endif()
    endif()

    include(ClangSanitizer)
endif() # COMPILER_IS_CLANG

### Intel compiler
if(COMPILER_IS_INTEL)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Set Intel release flags")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST")
endif() # COMPILER_IS_INTEL

# Profiling
if(OGS_PROFILE)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
        message(STATUS "When using profiling you should set CMAKE_BUILD_TYPE \
            to Release.")
    endif()
    set(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG")
    # clang compiler does not know the following flags
    if(NOT COMPILER_IS_CLANG)
        set(PROFILE_FLAGS "${PROFILE_FLAGS} -fno-inline-functions-called-once \
            -fno-optimize-sibling-calls")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
endif()

### Windows
if(WIN32)
    ## For Visual Studio compiler
    if(MSVC)
        if(OGS_CPU_ARCHITECTURE STREQUAL "native")
            set(CPU_FLAGS "/favor:blend")
        else()
            set(CPU_FLAGS "/favor:${OGS_CPU_ARCHITECTURE}")
        endif()
        if(OGS_ENABLE_AVX2)
            set(CPU_FLAGS "${CPU_FLAGS} /arch:AVX2")
        endif()
        add_definitions(
            -D_CRT_SECURE_NO_WARNINGS
            -D_CRT_NONSTDC_NO_WARNINGS
            -D_CRT_XNONSTDC_NO_WARNINGS
            -D__restrict__=__restrict   # this fixes #5
            # This fixes compile errors with
            # std::numeric_limits<T>::min() / max()
            -DNOMINMAX
        )
        # Sets warning level 3 and ignores some warnings
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPU_FLAGS} /W3 /wd4290 \
            /wd4267 /wd4996 /bigobj")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} \
            /ZI /Od /Ob0")

        DisableCompilerFlag(DEBUG /RTC1)
    # cygwin
    else()
        if(CMAKE_CROSSCOMPILING)
            message(STATUS "Crosscompiling for Windows with MinGW.")
        else()
            message(STATUS "Might be GCC under cygwin.")
        endif()
    endif()
endif()

# Missing OpenMP 3.0 implementation fix for Windows, this fixes #6
if(MSVC)
    add_definitions(-DOPENMP_LOOP_TYPE=int)
else()
    add_definitions(-DOPENMP_LOOP_TYPE=unsigned)
endif()

# Set additional user-given compiler flags. The given flags must follow the
# preceding cxx flags definition in order to override earlier flags, e.g. for
# optimization.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OGS_CXX_FLAGS}")

# Copy CXX to C flags. To have these flags in C libraries (e.g. metis) too.
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
