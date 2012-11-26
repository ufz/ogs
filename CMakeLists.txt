#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

# Project name
PROJECT( OGS-6 )

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 0)
SET (OGS_VERSION_MINOR 0)
SET (OGS_VERSION_PATCH 1)
SET (OGS_RELEASE_PERSONS "LB/TF/KR")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}")
SET (OGS_VERSION_AND_PERSONS "${OGS_VERSION} (${OGS_RELEASE_PERSONS})")
SET (OGS_DATE "2012-08-20")

### CMake includes ###
INCLUDE(scripts/cmake/CheckTypeSizes.cmake)
INCLUDE(scripts/cmake/FindIncludeHeader.cmake)
INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/SubmoduleSetup.cmake)
IF(NOT OGS_NO_EXTERNAL_LIBS)
	INCLUDE(scripts/cmake/Find.cmake)
ENDIF() # NOT OGS_NO_EXTERNAL_LIBS
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)
INCLUDE(scripts/cmake/Test.cmake)
IF(OGS_COVERAGE)
	INCLUDE(scripts/cmake/Coverage.cmake)
ENDIF()

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
OPTION(OGS_BUILD_UTILS "Should the utilities programms be built?" OFF)

OPTION(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

# Logging
OPTION(OGS_DISABLE_LOGGING "Disables all logog messages." OFF)

# Search paths
SET(OGS_LIBS_DIR "" CACHE PATH "The path to the compiled third party libs (mainly used for Windows).")

# Print CMake variable values
IF (OGS_CMAKE_DEBUG)
	INCLUDE(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
ENDIF (OGS_CMAKE_DEBUG)

OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
IF (OGS_BUILD_INFO)
	ADD_DEFINITIONS (-DOGS_BUILD_INFO)
ENDIF (OGS_BUILD_INFO)

# Code coverage
OPTION(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)

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
ADD_SUBDIRECTORY( ThirdParty )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/ThirdParty )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/ThirdParty/quickcheck )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/ThirdParty/zlib )

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
	ADD_SUBDIRECTORY( Tests )
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
IF(OGS_BUILD_UTILS)
	ADD_SUBDIRECTORY( Utils/SimpleMeshCreation )
ENDIF() # OGS_BUILD_UTILS

CONFIGURE_FILE (BaseLib/BuildInfo.h.in ${PROJECT_BINARY_DIR}/BaseLib/BuildInfo.h)
CONFIGURE_FILE (BaseLib/Configure.h.in ${PROJECT_BINARY_DIR}/BaseLib/Configure.h)
