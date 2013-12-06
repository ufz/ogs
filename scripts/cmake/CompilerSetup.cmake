INCLUDE(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
INCLUDE(SetDefaultBuildType)
INCLUDE(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Debug)
INCLUDE(MSVCMultipleProcessCompile) # /MP switch (multi processor) for VS

# Better Clang warning suppression, see http://www.openwalnut.org/issues/230
SET( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem" CACHE STRING "" FORCE )

# Set compiler helper variables
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(COMPILER_IS_CLANG TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(COMPILER_IS_GCC TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    SET(COMPILER_IS_INTEL TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(COMPILER_IS_MSVC TRUE)
ENDIF () # CMAKE_CXX_COMPILER_ID

# Set additional user-given compiler flags
SET(CMAKE_CXX_FLAGS ${OGS_CXX_FLAGS})

### GNU C/CXX compiler
IF(COMPILER_IS_GCC)
		get_gcc_version(GCC_VERSION)
		IF(GCC_VERSION VERSION_LESS "4.6")
			MESSAGE(FATAL_ERROR "GCC minimum required version is 4.6! You are using ${GCC_VERSION}.")
		ENDIF()
		IF(GCC_VERSION VERSION_LESS "4.7")
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
		ELSE()
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		ENDIF()
		IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
				MESSAGE(STATUS "Set GCC release flags")
				IF(APPLE AND GCC_VERSION VERSION_LESS "4.3" AND NOT "${CMAKE_GENERATOR}" STREQUAL "Xcode" )
					# -march=native does not work here when on normal gcc compiler
					# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=33144
					SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -msse4.2 -DNDEBUG")
				ELSE()
					SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -msse4.2 -DNDEBUG")
					# Disable -march=native on mac or Ninja generator
					IF(NOT APPLE AND NOT "${CMAKE_GENERATOR}" STREQUAL "Ninja")
						SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
					ENDIF()
				ENDIF()
		ENDIF()
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra")
ENDIF() # COMPILER_IS_GCC

IF(COMPILER_IS_CLANG)
	IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
		MESSAGE(FATAL_ERROR "Aborting: Clang 3.3 is required! Found version ${CMAKE_CXX_COMPILER_VERSION}")
	ENDIF()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Weverything -Wno-c++98-compat-pedantic")
ENDIF() # COMPILER_IS_CLANG

### Intel compiler
IF (COMPILER_IS_INTEL)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
				MESSAGE(STATUS "Set Intel release flags")
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
		ENDIF()
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -O3 -no-prec-div -DNDEBUG")
ENDIF() # COMPILER_IS_INTEL

# Profiling
IF (OGS_PROFILE)
	IF( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
		MESSAGE(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
	ENDIF()
	SET(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG")
	# clang compiler does not know the following flags
	IF(NOT COMPILER_IS_CLANG)
		SET(PROFILE_FLAGS "${PROFILE_FLAGS} -fno-inline-functions-called-once -fno-optimize-sibling-calls")
	ENDIF()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
ENDIF (OGS_PROFILE)

### Windows
IF (WIN32)
	## For Visual Studio compiler
	IF (MSVC)
		ADD_DEFINITIONS(
			-D_CRT_SECURE_NO_WARNINGS
			-D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS
			-D__restrict__=__restrict   # this fixes #5
			-DNOMINMAX # This fixes compile errors with std::numeric_limits<T>::min() / max()
		)
		# Sets warning level 3 and ignores some warnings
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")
		SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /ZI /Od /Ob0")

		DisableCompilerFlag(DEBUG /RTC1)
	# cygwin
	ELSE (MSVC)
		MESSAGE (STATUS "Might be GCC under cygwin.")
		ADD_DEFINITIONS( -DGCC )
	ENDIF (MSVC)
ENDIF (WIN32)

# Missing OpenMP 3.0 implementation fix for Windows, this fixes #6
IF(MSVC)
	ADD_DEFINITIONS(-DOPENMP_LOOP_TYPE=int)
ELSE()
	ADD_DEFINITIONS(-DOPENMP_LOOP_TYPE=unsigned)
ENDIF()
