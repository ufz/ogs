#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 3.1)

# Set CMake policies
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0054 NEW)

# Project name
project( OGS-6 )

###########################
### Preliminary Options ###
###########################
option(OGS_NO_BOOST_DOWNLOAD "Always use the system installed Boost?" OFF)
set(CMAKE_LIBRARY_SEARCH_PATH "" CACHE PATH
	"Additional library installation path, e.g. /opt/local or C:/libs")

### CMake includes ###
include(scripts/cmake/CheckTypeSizes.cmake)
include(scripts/cmake/Functions.cmake)
include(scripts/cmake/CMakeSetup.cmake)
include(scripts/cmake/CompilerSetup.cmake)
include(scripts/cmake/Find.cmake)
include(scripts/cmake/SubmoduleSetup.cmake)
include(scripts/cmake/ProjectSetup.cmake)
include(scripts/cmake/DocumentationSetup.cmake)
include(scripts/cmake/test/Test.cmake)
if(OGS_COVERAGE AND NOT IS_SUBPROJECT)
	include(scripts/cmake/Coverage.cmake)
endif()

####################
### More Options ###
####################

# Profiling
if((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	option(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
endif() # GCC AND GPROF_PATH

# Enable / Disable parts
option(OGS_DONT_USE_QT "Disables all Qt specific code." OFF)

option(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
option(OGS_BUILD_TESTS "Should the test executables be built?" ON)
option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
if(OGS_BUILD_GUI)
	add_definitions(-DOGS_BUILD_GUI)
endif()
option(OGS_BUILD_UTILS "Should the utilities programms be built?" OFF)

option(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

# Linear solvers
option(OGS_USE_LIS "Use Lis" OFF)

# Parallel computing: vector and matrix algebraic caculation, solvers
option(OGS_USE_PETSC "Use PETSc routines" OFF)

# Paralleization
option(OGS_USE_MPI "Use MPI" OFF)

# Eigen
option(OGS_USE_EIGEN "Use EIGEN for local matrix and vector" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" ON)

# Logging
option(OGS_DISABLE_LOGGING "Disables all logog messages." OFF)

# Compiler flags
set(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")

# Print CMake variable values
if (OGS_CMAKE_DEBUG)
	include(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
endif ()

# Code coverage
option(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)

# Packaging
option(OGS_DOWNLOAD_ADDITIONAL_CONTENT "Should addional content such as manuals be downloaded and packaged?" OFF)
include(scripts/cmake/packaging/Pack.cmake)

######################
### Subdirectories ###
######################

# External projects
option(Boost_USE_STATIC_LIBS "This option has to be set from the command line (before boost gets found)!" OFF)
include(scripts/cmake/ExternalProjectBoost.cmake)
option(OGS_SYSTEM_VTK_ONLY "If On VTK / ParaView will be searched if installed and will not be build locally." OFF)
include(scripts/cmake/ExternalProjectCatalyst.cmake)

include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
if(OGS_USE_EIGEN)
	include(scripts/cmake/ExternalProjectEigen.cmake)
endif()

# Add subdirectories with the projects
add_subdirectory( ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/gtest/include )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/autocheck/include )
include_directories( SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/zlib )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include )

if(OGS_USE_LIS)
    add_definitions(-DUSE_LIS)
    include_directories(SYSTEM ${LIS_INCLUDE_DIR})
endif()

# Use PETSc package for global assembly and linear sover for parallel computing
if(OGS_USE_PETSC)
    add_definitions(-DUSE_PETSC)
    set(OGS_USE_MPI ON)
endif()

# Use MPI
if(OGS_USE_MPI)
    add_definitions(-DUSE_MPI)
endif()

if(OGS_USE_EIGEN)
#	add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=std::size_t)
	add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO)
	add_definitions(-DOGS_USE_EIGEN)
	if (EIGEN_NO_DEBUG)
		add_definitions(-DEIGEN_NO_DEBUG)
	endif()
	include_directories (SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory( Applications )
add_subdirectory( AssemblerLib )
add_subdirectory( BaseLib )
# TODO This is a hack but we have to make sure that Boost is built first
if(TARGET Boost)
	add_dependencies(BaseLib Boost)
endif()
add_subdirectory( FileIO )
add_subdirectory( GeoLib )
add_subdirectory( InSituLib )
add_subdirectory( MathLib )
add_subdirectory( MeshLib )
add_subdirectory( MeshGeoToolsLib )
add_subdirectory( NumLib )
add_subdirectory( ProcessLib )
if( OGS_BUILD_TESTS AND NOT IS_SUBPROJECT )
	add_subdirectory( Tests )

	if(OGS_USE_MPI)
		add_subdirectory( SimpleTests/MeshTests/MPI )
	else()
		add_subdirectory( SimpleTests/MatrixTests EXCLUDE_FROM_ALL )
		add_subdirectory( SimpleTests/MeshTests EXCLUDE_FROM_ALL )
		if(NOT MSVC AND BLAS_FOUND AND LAPACK_FOUND)
			add_subdirectory( SimpleTests/SolverTests EXCLUDE_FROM_ALL )
		endif()
	endif()

	# Create a target 'data', which downloads all referenced data sets into the build tree
	# This has to be defined after all tests are defined
	ExternalData_Add_Target(data)
	set_target_properties(data PROPERTIES EXCLUDE_FROM_ALL 1)
	add_dependencies(ctest data)
endif() # OGS_BUILD_TESTS

# The configuration must be called from the source dir and not BaseLib/.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/BaseLib/BuildInfo.cpp.in"
	"${CMAKE_CURRENT_BINARY_DIR}/BaseLib/BuildInfo.cpp" @ONLY)
