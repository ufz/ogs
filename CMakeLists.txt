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

list(APPEND CMAKE_PREFIX_PATH
    $ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
    $ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
    ${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

###########################
### Preliminary Options ###
###########################
set(CMAKE_LIBRARY_SEARCH_PATH "" CACHE PATH
    "Additional library installation path, e.g. /opt/local or C:/libs")
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, defaults to native.")
option(OGS_BUILD_TESTS "Should the test executables be built?" ON)

### CMake includes ###
include(scripts/cmake/ConanSetup.cmake)
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
option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
if(OGS_BUILD_GUI)
    add_definitions(-DOGS_BUILD_GUI)
endif()
option(OGS_BUILD_UTILS "Should the utilities programms be built?" OFF)

if(WIN32)
    option(OGS_BUILD_SWMM "Should the SWMM interface be built?" OFF)
endif()
if(NOT WIN32 AND OGS_BUILD_SWMM)
    message(FATAL_ERROR "OGS_BUILD_SWMM requires Windows!")
endif()

option(OGS_BUILD_METIS "Should metis and the partmesh util be built?" OFF)

option(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

# Linear solvers
option(OGS_USE_LIS "Use Lis" OFF)

# Parallel computing: vector and matrix algebraic caculation, solvers
option(OGS_USE_PETSC "Use PETSc routines" OFF)

# Paralleization
option(OGS_USE_MPI "Use MPI" OFF)

# Eigen
option(OGS_USE_EIGEN "Use Eigen linear solver" ON)
option(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Use dynamically allocated shape matrices" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# Logging
option(OGS_DISABLE_LOGGING "Disables all logog messages." OFF)

# Debug
option(OGS_FATAL_ABORT "Abort in OGS_FATAL" OFF)

# Compiler flags
set(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")
option(STL_NO_DEBUG "Disable STL debug in debug build" OFF)

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

# Third-party libraries
# Default checks for system first, then builds locally
set(OGS_LIB_VTK "Default" CACHE STRING "Which VTK library should be used?")
set(OGS_LIB_BOOST "Default" CACHE STRING "Which Boost library should be used?")
set(OGS_LIB_EIGEN "Default" CACHE STRING "Which Eigen library should be used?")
set_property(CACHE
    OGS_LIB_VTK
    OGS_LIB_BOOST
    OGS_LIB_EIGEN
    PROPERTY STRINGS "Default" "System" "Local")

# Options controlling which FEM elements will be compiled
set(OGS_MAX_ELEMENT_DIM   3 CACHE STRING "Maximum dimension of FEM elements to be built.")
set(OGS_MAX_ELEMENT_ORDER 2 CACHE STRING "Maximum order of FEM elements to be built.")
option(OGS_ENABLE_ELEMENT_SIMPLEX "Build FEM elements for simplices (triangles, tetrahedra)." ON)
option(OGS_ENABLE_ELEMENT_CUBOID  "Build FEM elements for cuboids (quads, hexahedra)." ON)
option(OGS_ENABLE_ELEMENT_PRISM   "Build FEM elements for prisms." ON)
option(OGS_ENABLE_ELEMENT_PYRAMID "Build FEM elements for pyramids." ON)

###################
### Definitions ###
###################
if(OGS_USE_LIS)
    add_definitions(-DUSE_LIS)
    include_directories(SYSTEM ${LIS_INCLUDE_DIR})
    set(OGS_USE_EIGEN ON)
endif()

if(OGS_USE_PETSC)
    add_definitions(-DUSE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()

# Use MPI
if(OGS_USE_MPI)
    add_definitions(-DUSE_MPI)
endif()

add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO) # TODO check if needed
if (EIGEN_NO_DEBUG)
    add_definitions(-DEIGEN_NO_DEBUG)
endif()
if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
    add_definitions(-DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
endif()

if(OGS_USE_EIGEN)
    add_definitions(-DOGS_USE_EIGEN)
endif()

if (OGS_FATAL_ABORT)
    add_definitions(-DOGS_FATAL_ABORT)
endif()

if(OGS_BUILD_TESTS)
    set(Data_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
    set(Data_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")
endif()

######################
### Subdirectories ###
######################

# External projects
include(scripts/cmake/ExternalProjectBoost.cmake)
include(VtkModules)
if(OGS_INSITU)
    include(ExternalProjectCatalyst)
else()
    include(ExternalProjectVtk)
endif()

include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
include(scripts/cmake/ExternalProjectEigen.cmake)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Add subdirectories with the projects
add_subdirectory( ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/gtest/include )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/autocheck/include )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include )
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/vtkGUISupportQt)

add_subdirectory( Applications )
add_subdirectory( BaseLib )
# TODO This is a hack but we have to make sure that Boost is built first
if(TARGET Boost)
    add_dependencies(BaseLib Boost)
endif()
add_subdirectory( GeoLib )
add_subdirectory( MaterialLib )
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
        add_subdirectory( SimpleTests/MatrixTests )
        add_subdirectory( SimpleTests/MeshTests )
    endif()
endif() # OGS_BUILD_TESTS

# The configuration must be called from the source dir and not BaseLib/.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/BaseLib/BuildInfo.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BaseLib/BuildInfo.cpp" @ONLY)

include(scripts/cmake/MarkVariablesAdvanced.cmake)
include(scripts/cmake/packaging/GetDependencies.cmake)
