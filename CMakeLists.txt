#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

# Project name
PROJECT( OGS-6 )

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 6)
SET (OGS_VERSION_MINOR 0)
SET (OGS_VERSION_PATCH "0 Beta")
SET (OGS_RELEASE_PERSONS "Core developers")
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
INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/Find.cmake)
INCLUDE(scripts/cmake/SubmoduleSetup.cmake)
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)
IF(CMAKE_VERSION VERSION_GREATER 2.8.10)
	INCLUDE(scripts/cmake/test/Test.cmake)
ENDIF()
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
OPTION(OGS_BUILD_VIS "Should in-situ visulization be built?" OFF)

OPTION(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

# Linear solvers
OPTION(OGS_USE_LIS "Use Lis" OFF)
OPTION(OGS_USE_PETSC "Use PETSc routines" OFF)

# Paralleization
OPTION(OGS_USE_MPI "Use MPI" OFF)

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
OPTION(OGS_SYSTEM_VTK_ONLY "If On VTK / ParaView will be searched if installed and will not be build locally." OFF)
INCLUDE(scripts/cmake/ExternalProjectCatalyst.cmake)

INCLUDE_DIRECTORIES( SYSTEM ${Boost_INCLUDE_DIRS} )
IF(OGS_USE_EIGEN)
	INCLUDE(scripts/cmake/ExternalProjectEigen.cmake)
ENDIF()

# Add subdirectories with the projects
ADD_SUBDIRECTORY( ThirdParty )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/gtest/include )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/quickcheck )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/zlib )

IF(OGS_USE_LIS)
    ADD_DEFINITIONS(-DUSE_LIS)
    INCLUDE_DIRECTORIES(SYSTEM ${LIS_INCLUDE_DIR})
ENDIF()

IF(OGS_USE_PETSC)
    ADD_DEFINITIONS(-DUSE_PETSC)
    SET(OGS_USE_MPI ON)
ENDIF()

IF(OGS_USE_MPI)
    ADD_DEFINITIONS(-DUSE_MPI)
ENDIF()

IF(OGS_USE_EIGEN)
#	ADD_DEFINITIONS(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=std::size_t)
	ADD_DEFINITIONS(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO)
	ADD_DEFINITIONS(-DOGS_USE_EIGEN)
	IF (EIGEN_NO_DEBUG)
		ADD_DEFINITIONS(-DEIGEN_NO_DEBUG)
	ENDIF()
	INCLUDE_DIRECTORIES (SYSTEM ${EIGEN3_INCLUDE_DIR})
ENDIF()

IF(OGS_BUILD_VIS)
	ADD_SUBDIRECTORY(InSituLib)
ENDIF()

ADD_SUBDIRECTORY( Applications )
ADD_SUBDIRECTORY( AssemblerLib )
ADD_SUBDIRECTORY( BaseLib )
# TODO This is a hack but we have to make sure that Boost is built first
IF(TARGET Boost)
	ADD_DEPENDENCIES(BaseLib Boost)
ENDIF()
ADD_SUBDIRECTORY( FileIO )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( InSituLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( MeshLib )
ADD_SUBDIRECTORY( MeshGeoToolsLib )
ADD_SUBDIRECTORY( NumLib )
ADD_SUBDIRECTORY( ProcessLib )
IF( OGS_BUILD_TESTS AND NOT IS_SUBPROJECT )
	IF(CMAKE_VERSION VERSION_LESS 2.8.11)
		MESSAGE(FATAL_ERROR "CMAKE 2.8.11 or higher is required to build the tests!")
	ENDIF()
	ADD_SUBDIRECTORY( Tests )
	ADD_SUBDIRECTORY( SimpleTests/MatrixTests EXCLUDE_FROM_ALL )
	ADD_SUBDIRECTORY( SimpleTests/MeshTests EXCLUDE_FROM_ALL )
	IF(NOT MSVC AND BLAS_FOUND AND LAPACK_FOUND)
		ADD_SUBDIRECTORY( SimpleTests/SolverTests EXCLUDE_FROM_ALL )
	ENDIF()
	# Create a target 'data', which downloads all referenced data sets into the build tree
	# This has to be defined after all tests are defined
	ExternalData_Add_Target(data)
	SET_TARGET_PROPERTIES(data PROPERTIES EXCLUDE_FROM_ALL 1)
	ADD_DEPENDENCIES(ctest data)
ENDIF() # OGS_BUILD_TESTS

# The configuration must be called from the source dir and not BaseLib/.
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/BaseLib/BuildInfo.cpp.in"
	"${CMAKE_CURRENT_BINARY_DIR}/BaseLib/BuildInfo.cpp" @ONLY)
