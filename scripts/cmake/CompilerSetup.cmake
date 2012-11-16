INCLUDE(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
INCLUDE(SetDefaultBuildType)
INCLUDE(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Debug)
INCLUDE(MSVCMultipleProcessCompile) # /MP switch (multi processor) for VS

# Set compiler helper variables
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	SET(COMPILER_IS_CLANG 1)
ENDIF ()
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET(COMPILER_IS_GCC 1)
ENDIF ()
IF (${CMAKE_C_COMPILER} MATCHES "icc.*$" OR ${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
	SET(COMPILER_IS_INTEL)
ENDIF ()

### GNU C/CXX compiler
IF(COMPILER_IS_GCC)
		get_gcc_version(GCC_VERSION)
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
					SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native -msse4.2 -DNDEBUG")
				ENDIF()
		ENDIF()
		# -g
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra")
		IF(COMPILER_IS_CLANG)
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
		ELSE()
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-nonansi-builtins")
		ENDIF() # COMPILER_IS_CLANG
		ADD_DEFINITIONS( -DGCC -Wfatal-errors)
ENDIF() # COMPILER_IS_GCC

### Intel compiler
IF (COMPILER_IS_INTEL)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
				MESSAGE(STATUS "Set Intel release flags")
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
		ENDIF()
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -O3 -no-prec-div -static -DNDEBUG")
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
