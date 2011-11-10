# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)

# Project name
project( OGS-6 )

# Set cmake module path 
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


# Set cmake module path 
#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeConfiguration")

### For GNU C/CXX
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set GCC release flags")
                SET(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -msse4.2 -DNDEBUG")
        ENDIF()
        # -g
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        ADD_DEFINITIONS(
                -DGCC
        )
ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

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

