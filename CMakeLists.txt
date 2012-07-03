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
IF(NOT OGS_NO_EXTERNAL_LIBS)
	INCLUDE(scripts/cmake/Find.cmake)
ENDIF() # NOT OGS_NO_EXTERNAL_LIBS
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)

###############
### Options ###
###############

# Profiling
IF((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
ENDIF() # GCC AND GPROF_PATH

# Disable Qt
OPTION(OGS_DONT_USE_QT "Disables all Qt specific code." OFF)

# Print CMake variable values
IF (OGS_CMAKE_DEBUG)
	INCLUDE(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
ENDIF (OGS_CMAKE_DEBUG)

OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
IF (OGS_BUILD_INFO)
	ADD_DEFINITIONS (-DOGS_BUILD_INFO)
ENDIF (OGS_BUILD_INFO)

OPTION(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

######################
### Subdirectories ###
######################

# Add subdirectories with the projects
ADD_SUBDIRECTORY( BaseLib )
ADD_SUBDIRECTORY( FemLib )
ADD_SUBDIRECTORY( FileIO )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( MeshLib )
ADD_SUBDIRECTORY( OGS )
ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
ADD_SUBDIRECTORY( SimpleTests/MeshTests )
IF(NOT MSVC)
        ADD_SUBDIRECTORY( SimpleTests/SolverTests )
ENDIF(NOT MSVC)

CONFIGURE_FILE (BaseLib/BuildInfo.h.in ${PROJECT_BINARY_DIR}/BaseLib/BuildInfo.h)
