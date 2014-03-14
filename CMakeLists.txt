#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

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


###########################
### Preliminary Options ###
###########################
OPTION(OGS_NO_BOOST_DOWNLOAD "Always use the system installed Boost?" OFF)
SET(CMAKE_LIBRARY_SEARCH_PATH "" CACHE PATH
	"Additional library installation path, e.g. /opt/local or C:/libs")

### CMake includes ###
INCLUDE(scripts/cmake/CheckTypeSizes.cmake)
INCLUDE(scripts/cmake/FindIncludeHeader.cmake)
INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/Find.cmake)
INCLUDE(scripts/cmake/SubmoduleSetup.cmake)
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)
INCLUDE(scripts/cmake/Test.cmake)
IF(OGS_COVERAGE AND NOT IS_SUBPROJECT)
	INCLUDE(scripts/cmake/Coverage.cmake)
ENDIF()

####################
### More Options ###
####################

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

# Linear solvers
OPTION(OGS_USE_LIS "Use Lis" OFF)
OPTION(OGS_USE_PETSC "Use PETSc routines" OFF)

# Eigen
OPTION(OGS_USE_EIGEN "Use EIGEN for local matrix and vector" ON)
OPTION(EIGEN_NO_DEBUG "Disables Eigen's assertions" ON)

# Logging
OPTION(OGS_DISABLE_LOGGING "Disables all logog messages." OFF)

# Compiler flags
SET(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")

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
IF(WIN32)
	OPTION(OGS_PACKAGING_NSIS "Do you want to package with NSIS?")
ENDIF()
INCLUDE(scripts/cmake/Pack.cmake)

######################
### Subdirectories ###
######################

# External projects
OPTION(Boost_USE_STATIC_LIBS "This option has to be set from the command line (before boost gets found)!" OFF)
INCLUDE(scripts/cmake/ExternalProjectBoost.cmake)
OPTION(OGS_SYSTEM_VTK_ONLY "If On VTK will be searched if installed and will not be build locally." OFF)
INCLUDE(scripts/cmake/ExternalProjectVtk.cmake)
INCLUDE_DIRECTORIES( SYSTEM ${Boost_INCLUDE_DIRS} )
IF(OGS_USE_EIGEN)
	INCLUDE(scripts/cmake/ExternalProjectEigen.cmake)
ENDIF()

# Add subdirectories with the projects
ADD_SUBDIRECTORY( ThirdParty )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/quickcheck )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/zlib )

IF(OGS_USE_LIS)
    ADD_DEFINITIONS(-DUSE_LIS)
ENDIF()

IF(OGS_USE_PETSC)
    ADD_DEFINITIONS(-DUSE_PETSC)
ENDIF()

IF(OGS_USE_EIGEN)
#	ADD_DEFINITIONS(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=std::size_t)
	ADD_DEFINITIONS(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO)
	ADD_DEFINITIONS(-DOGS_USE_EIGEN)
	IF (EIGEN_NO_DEBUG)
		ADD_DEFINITIONS(-DEIGEN_NO_DEBUG)
	ENDIF()
ENDIF()

IF(OGS_BUILD_GUI)
	ADD_DEFINITIONS(-DOGS_BUILD_GUI)
	ADD_SUBDIRECTORY(Gui)
ENDIF() # OGS_BUILD_GUI

ADD_SUBDIRECTORY( AssemblerLib )
ADD_SUBDIRECTORY( BaseLib )
# TODO This is a hack but we have to make sure that Boost is built first
ADD_DEPENDENCIES(BaseLib Boost)
ADD_SUBDIRECTORY( FileIO )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( MeshLib )
ADD_SUBDIRECTORY( MeshGeoToolsLib )
ADD_SUBDIRECTORY( NumLib )
IF( OGS_BUILD_CLI )
	ADD_SUBDIRECTORY( OGS )
ENDIF() # OGS_BUILD_CLI
IF( OGS_BUILD_TESTS AND NOT IS_SUBPROJECT )
	ADD_SUBDIRECTORY( Tests )
	ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
	ADD_SUBDIRECTORY( SimpleTests/MeshTests )
	IF(NOT MSVC AND BLAS_FOUND AND LAPACK_FOUND)
		ADD_SUBDIRECTORY( SimpleTests/SolverTests )
	ENDIF()
ENDIF() # OGS_BUILD_TESTS
IF(OGS_BUILD_UTILS AND NOT IS_SUBPROJECT)
	ADD_SUBDIRECTORY( Utils/SimpleMeshCreation )
	ADD_SUBDIRECTORY( Utils/FileConverter )
	ADD_SUBDIRECTORY( Utils/MeshEdit )
ENDIF() # OGS_BUILD_UTILS

CONFIGURE_FILE (BaseLib/BuildInfo.h.in ${PROJECT_BINARY_DIR}/BaseLib/BuildInfo.h)
CONFIGURE_FILE (BaseLib/Configure.h.in ${PROJECT_BINARY_DIR}/BaseLib/Configure.h)
