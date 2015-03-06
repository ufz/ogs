include(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
include(SetDefaultBuildType)
include(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Debug)
include(MSVCMultipleProcessCompile) # /MP switch (multi processor) for VS
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Set compiler helper variables

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(COMPILER_IS_CLANG TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(COMPILER_IS_GCC TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(COMPILER_IS_INTEL TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    set(COMPILER_IS_MSVC TRUE)
endif() # CMAKE_CXX_COMPILER_ID

# Better Clang warning suppression, see http://www.openwalnut.org/issues/230
if(NOT COMPILER_IS_MSVC)
	set( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem" CACHE STRING "" FORCE )
endif()

# Set additional user-given compiler flags
set(CMAKE_CXX_FLAGS ${OGS_CXX_FLAGS})

### GNU C/CXX compiler
if(COMPILER_IS_GCC)
	get_gcc_version(GCC_VERSION)
	if(GCC_VERSION VERSION_LESS "4.8")
		message(FATAL_ERROR "GCC minimum required version is 4.8! You are using ${GCC_VERSION}.")
	endif()
	if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
		message(STATUS "Set GCC release flags")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -msse4.2 -DNDEBUG")
		# Disable -march=native on mac or Ninja generator
		if(NOT APPLE AND NOT "${CMAKE_GENERATOR}" STREQUAL "Ninja")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
		endif()
	else()
		# Enable assertions in STL in debug mode.
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_ASSERT -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -Wall -Wextra")
endif() # COMPILER_IS_GCC

### Clang
if(COMPILER_IS_CLANG)
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
		message(FATAL_ERROR "Aborting: Clang 3.3 is required! Found version ${CMAKE_CXX_COMPILER_VERSION}")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-c++98-compat-pedantic")
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		# Enable assertions in STL in debug mode.
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_ASSERT -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY")
	endif()
endif() # COMPILER_IS_CLANG

### Intel compiler
if(COMPILER_IS_INTEL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
		message(STATUS "Set Intel release flags")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -O3 -no-prec-div -DNDEBUG")
endif() # COMPILER_IS_INTEL

# Profiling
if(OGS_PROFILE)
	if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
		message(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
	endif()
	set(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG")
	# clang compiler does not know the following flags
	if(NOT COMPILER_IS_CLANG)
		set(PROFILE_FLAGS "${PROFILE_FLAGS} -fno-inline-functions-called-once -fno-optimize-sibling-calls")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
endif ()

### Windows
if (WIN32)
	## For Visual Studio compiler
	if (MSVC)
		add_definitions(
			-D_CRT_SECURE_NO_WARNINGS
			-D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS
			-D__restrict__=__restrict   # this fixes #5
			-DNOMINMAX # This fixes compile errors with std::numeric_limits<T>::min() / max()
		)
		# Sets warning level 3 and ignores some warnings
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267 /wd4996")
		set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /ZI /Od /Ob0")

		DisableCompilerFlag(DEBUG /RTC1)
	# cygwin
	else()
		message (STATUS "Might be GCC under cygwin.")
		add_definitions( -DGCC )
	endif()
endif ()

# Missing OpenMP 3.0 implementation fix for Windows, this fixes #6
if(MSVC)
	add_definitions(-DOPENMP_LOOP_TYPE=int)
else()
	add_definitions(-DOPENMP_LOOP_TYPE=unsigned)
endif()
