#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)

# Project name
project( OGS-6 )

### CMake includes ###
INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/Find.cmake)
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)

###############
### Options ###
###############

# Profiling
IF((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
ENDIF() # GCC AND GPROF_PATH

######################
### Subdirectories ###
######################

# Add subdirectories with the projects
ADD_SUBDIRECTORY( Base )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( MeshLib )
ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
IF(NOT MSVC)
        ADD_SUBDIRECTORY( SimpleTests/SolverTests )
ENDIF(NOT MSVC)

