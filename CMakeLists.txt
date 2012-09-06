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

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 0)
SET (OGS_VERSION_MINOR 0)
SET (OGS_VERSION_PATCH 1)
SET (OGS_RELEASE_PERSONS "LB/TF/KR")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
SET (OGS_DATE "2012-08-20")

###############
### Options ###
###############

# Profiling
IF((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
ENDIF() # GCC AND GPROF_PATH

# Enable / Disable parts
OPTION(OGS_DONT_USE_QT "Disables all Qt specific code." OFF)

OPTION(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
OPTION(OGS_BUILD_TESTS "Should the test executables be built?" ON)
OPTION(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)

OPTION(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

# Print CMake variable values
IF (OGS_CMAKE_DEBUG)
	INCLUDE(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
ENDIF (OGS_CMAKE_DEBUG)

OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
IF (OGS_BUILD_INFO)
	ADD_DEFINITIONS (-DOGS_BUILD_INFO)
ENDIF (OGS_BUILD_INFO)

# Packaging
OPTION(OGS_PACKAGING "Creating installers / packages" OFF)
OPTION_REQUIRES(OGS_PACKAGING_ZIP "Do you want to package as zip?" OGS_PACKAGING)
IF(OGS_PACKAGING)
	INCLUDE( scripts/cmake/Pack.cmake)
ENDIF() #OGS_PACKAGING

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
IF( OGS_BUILD_CLI )
	ADD_SUBDIRECTORY( OGS )
ENDIF() # OGS_BUILD_CLI
IF( OGS_BUILD_TESTS )
	ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
	ADD_SUBDIRECTORY( SimpleTests/MeshTests )
	IF(NOT MSVC AND BLAS_FOUND AND LAPACK_FOUND)
		ADD_SUBDIRECTORY( SimpleTests/SolverTests )
	ENDIF()
ENDIF() # OGS_BUILD_TESTS
IF(OGS_BUILD_GUI)
	ADD_DEFINITIONS(-DOGS_BUILD_GUI)
	ADD_SUBDIRECTORY(Gui)
ENDIF() # OGS_BUILD_GUI

CONFIGURE_FILE (BaseLib/BuildInfo.h.in ${PROJECT_BINARY_DIR}/BaseLib/BuildInfo.h)
CONFIGURE_FILE (BaseLib/Configure.h.in ${PROJECT_BINARY_DIR}/BaseLib/Configure.h)
