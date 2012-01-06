# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)

# Project name
project( OGS-6 )

### CMake setup ###
INCLUDE(scripts/cmake/CMakeSetup.cmake)

# Finds
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

### For GNU C/CXX
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
		get_gcc_version(GCC_VERSION)
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set GCC release flags")
				IF(APPLE AND GCC_VERSION VERSION_LESS "4.3" AND NOT "${CMAKE_GENERATOR}" STREQUAL "Xcode" )
					# -march=native does not work here when on normal gcc compiler
					# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=33144
					SET(CMAKE_CXX_FLAGS "-O3 -mtune=native -msse4.2 -DNDEBUG")
				ELSE()
                	SET(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -msse4.2 -DNDEBUG")
				ENDIF()
        ENDIF()
        # -g
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")
        ADD_DEFINITIONS(
                -DGCC
        )
ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

IF (${CMAKE_C_COMPILER} MATCHES "icc.*$" OR ${CMAKE_CXX_COMPILER} MATCHES "icpc.*$") 
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set Intel release flags")
                SET(CMAKE_CXX_FLAGS "-O3 -DNDEBUG")
        ENDIF()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
ENDIF(${CMAKE_C_COMPILER} MATCHES "icc.*$" OR ${CMAKE_CXX_COMPILER} MATCHES
"icpc.*$")


# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

# Add subdirectories with the projects
ADD_SUBDIRECTORY( Base )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
IF(NOT MSVC)
        ADD_SUBDIRECTORY( SimpleTests/SolverTests )
ENDIF(NOT MSVC)

