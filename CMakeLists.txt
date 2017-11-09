#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 3.1)

# Set CMake policies
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0054 NEW)

if(WIN32 AND
    NOT ${CMAKE_GENERATOR} STREQUAL Ninja AND
    NOT ${CMAKE_GENERATOR} MATCHES 2015) # does not work on 2015
    # Force 64-bit VS compiler
    set(CMAKE_GENERATOR_TOOLSET host=x64)
endif()

# Project name
project( OGS-6 )

list(APPEND CMAKE_PREFIX_PATH
    $ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
    $ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
    ${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

###########################
### Preliminary Options ###
###########################
option(BUILD_SHARED_LIBS "Create shared libraries?" OFF)
set(CMAKE_LIBRARY_SEARCH_PATH "" CACHE PATH
    "Additional library installation path, e.g. /opt/local or C:/libs")
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, \
    defaults to native (*nix) / blend (MSVC).")
option(OGS_ENABLE_AVX2 "Enable the use of AVX2 instructions" OFF)
option(OGS_BUILD_TESTS "Should the test executables be built?" ON)
option(OGS_USE_PCH "Should pre-compiled headers be used?" ON)
option(OGS_USE_CONAN "Should Conan package manager be used?" OFF)
set(OGS_CONAN_BUILD "missing" CACHE STRING "Possible values: all, missing, \
    never or list of libs to build" )

### CMake includes ###
include(scripts/cmake/PreFind.cmake)
include(scripts/cmake/CheckTypeSizes.cmake)
include(scripts/cmake/Functions.cmake)
include(scripts/cmake/SubmoduleSetup.cmake)
include(scripts/cmake/CMakeSetup.cmake)
include(scripts/cmake/ConanSetup.cmake)
include(scripts/cmake/CompilerSetup.cmake)
include(scripts/cmake/Find.cmake)
include(scripts/cmake/CCacheSetup.cmake)
include(scripts/cmake/ProjectSetup.cmake)
include(scripts/cmake/DocumentationSetup.cmake)
include(scripts/cmake/test/Test.cmake)
if(OGS_COVERAGE AND NOT IS_SUBPROJECT)
    include(scripts/cmake/Coverage.cmake)
endif()
include(scripts/cmake/Web.cmake)

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

option(OGS_INSITU "Builds OGS with insitu visualization capabilities." OFF)

# Linear solvers
option(OGS_USE_LIS "Use Lis" OFF)

# Parallel computing: vector and matrix algebraic caculation, solvers
option(OGS_USE_PETSC "Use PETSc routines" OFF)

# Eigen
option(OGS_USE_EIGEN "Use Eigen linear solver" ON)
option(OGS_USE_EIGEN_UNSUPPORTED "Use Eigen unsupported modules" ON)
option(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN "" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# We assume that it's save to use vectorization with Eigen < 3.3 (strictly smaller than 3.3.!).
# At least we only observed vectorization issues with Eigen 3.3.x.
# If you want to use Eigen vectorization, make sure that you run all the ctests several times, e.g.:
# $ ctest --repeat-until-fail 50
# You might also want to take a look at https://github.com/ufz/ogs/issues/1881.
option(EIGEN_DONT_VECTORIZE "Disables explicit vectorization when defined." ON)

set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Default" CACHE STRING "Use dynamically allocated shape matrices")
set_property(CACHE OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
             PROPERTY STRINGS "Default" "ON" "OFF")

if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES STREQUAL "Default")
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR
        CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL OFF)
    else()
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL ON)
    endif()
else()
    set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL ${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES})
endif()

# MKL
option(OGS_USE_MKL "Use Intel MKL" OFF)

# Logging
set(OGS_LOG_LEVEL "LOGOG_LEVEL_DEBUG" CACHE STRING "Set logging level included in compilation.")
set_property(CACHE OGS_LOG_LEVEL PROPERTY STRINGS
    LOGOG_LEVEL_NONE
    LOGOG_LEVEL_ERROR
    LOGOG_LEVEL_WARN
    LOGOG_LEVEL_INFO
    LOGOG_LEVEL_DEBUG
    LOGOG_LEVEL_ALL
)

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

option(OGS_SYNC_SUBMODULES "Sync git submodules on CMake run. Overwrites custom URLs!" ON)

# Packaging
option(OGS_DOWNLOAD_ADDITIONAL_CONTENT "Should additional content such as manuals be downloaded and packaged?" OFF)
option(OGS_PACKAGE_DEPENDENCIES "Packages shared libraries when enabled." OFF)
include(scripts/cmake/packaging/Pack.cmake)

# Third-party libraries, names come from Conan package names
set(OGS_LIBS VTK Boost Eigen3 Qt Shapelib libgeotiff CACHE INTERNAL "")
foreach(LIB ${OGS_LIBS})
    set(OGS_LIB_${LIB} "Default" CACHE STRING "Which ${LIB} library should be used?")
    set_property(CACHE OGS_LIB_${LIB} PROPERTY STRINGS "Default" "System")
endforeach(LIB)

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

if(OGS_USE_MKL)
    add_definitions(-DUSE_MKL)
    include_directories(SYSTEM ${MKL_INCLUDE_DIR})
endif()

if(OGS_USE_PETSC)
    add_definitions(-DUSE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()

# Use MPI
if(OGS_USE_MPI)
    add_definitions(-DUSE_MPI)
endif()

# Eigen
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()

if (EIGEN_NO_DEBUG)
    add_definitions(-DEIGEN_NO_DEBUG)
endif()

if (EIGEN_DONT_VECTORIZE)
    add_definitions(-DEIGEN_DONT_VECTORIZE)
endif()

if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL)
    add_definitions(-DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
endif()

if(OGS_USE_EIGEN)
    add_definitions(-DOGS_USE_EIGEN)
    if(OGS_USE_EIGEN_UNSUPPORTED)
        add_definitions(-DUSE_EIGEN_UNSUPPORTED)
    endif()
    if(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN)
        add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
    endif()
endif()

if(MSVC AND OGS_32_BIT)
    add_definitions(-DEIGEN_MAX_ALIGN_BYTES=0 -DEIGEN_DONT_ALIGN)
endif()
# End Eigen

if (OGS_FATAL_ABORT)
    add_definitions(-DOGS_FATAL_ABORT)
endif()

if(OGS_BUILD_TESTS)
    set(Data_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Tests/lfs-data CACHE INTERNAL "")
    set(Data_BINARY_DIR ${PROJECT_BINARY_DIR}/Tests/lfs-data CACHE INTERNAL "")
endif()

# Logging level
add_definitions(-DLOGOG_LEVEL=${OGS_LOG_LEVEL})


######################
### Subdirectories ###
######################

# External projects
if(OGS_INSITU)
    include("${PARAVIEW_USE_FILE}")
    add_definitions(-DUSE_INSITU)
endif()


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
add_subdirectory( MathLib )
add_subdirectory( MeshLib )
add_subdirectory( MeshGeoToolsLib )
add_subdirectory( NumLib )
if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR OGS_BUILD_TESTS)
    add_subdirectory( MaterialLib )
    add_subdirectory( ProcessLib )
endif()
if( OGS_BUILD_TESTS AND NOT IS_SUBPROJECT )
    add_subdirectory( Tests )

    if(OGS_USE_MPI)
        add_subdirectory( SimpleTests/MeshTests/MPI )
    else()
        add_subdirectory( SimpleTests/MeshTests )
    endif()
endif() # OGS_BUILD_TESTS

# The configuration must be called from the source dir and not BaseLib/.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/BaseLib/BuildInfo.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BaseLib/BuildInfo.cpp" @ONLY)

include(scripts/cmake/MarkVariablesAdvanced.cmake)
unset(PRE_INSTALL_RUN CACHE)
