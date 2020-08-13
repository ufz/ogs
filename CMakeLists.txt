# ---- OGS-6 Project ----
cmake_minimum_required(VERSION 3.14)

# Set CMake policies
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.15)
    cmake_policy(SET CMP0092 NEW)
endif()

project(OGS-6)

include(scripts/cmake/CMakeSetup.cmake)
include(ParseCMakeArgs)
include(CTest)

# ---- Preliminary Options ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "Create shared libraries?" OFF)
option(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
set(CMAKE_LIBRARY_SEARCH_PATH
    ""
    CACHE PATH
          "Additional library installation path, e.g. /opt/local or C:/libs")
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, \
    defaults to native (*nix) / blend (MSVC).")
option(OGS_ENABLE_AVX2 "Enable the use of AVX2 instructions" OFF)
option(OGS_USE_CONAN "Should Conan package manager be used?" ON)
set(OGS_CONAN_BUILD "missing" CACHE STRING "Possible values: all, missing, \
    never or list of libs to build")
if(WIN32)
    option(OGS_DISABLE_CLCACHE "Disables clcache compiler cache." OFF)
else()
    option(OGS_DISABLE_CCACHE "Disables ccache compiler cache." OFF)
endif()
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.16)
    option(OGS_USE_UNITY_BUILDS "Enables Unity builds for faster compilation." ON)
endif()

# Third-party libraries, names come from Conan package names
set(OGS_LIBS
    vtk
    boost
    eigen
    qt
    shapelib
    libgeotiff
    petsc
    lis
    cvode
    tfel
    CACHE INTERNAL "")
foreach(LIB ${OGS_LIBS})
    set(OGS_LIB_${LIB}
        "Default"
        CACHE STRING "Which ${LIB} library should be used?")
    set_property(CACHE OGS_LIB_${LIB} PROPERTY STRINGS "Default" "System")
endforeach(LIB)

if(OGS_USE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()
option(OGS_USE_CVODE "Use the Sundials CVODE module?" OFF)

# ---- CMake includes ----
include(PreFind)
include(SubmoduleSetup)
include(ProcessesSetup)
include(ProjectSetup)
include(Versions)
include(Functions)
include(ConanSetup)
include(CompilerSetup)
include(Find)
include(CCacheSetup)
include(CLCacheSetup)
include(DocumentationSetup)
include(test/Test)
if(OGS_COVERAGE AND NOT IS_SUBPROJECT)
    include(Coverage)
endif()
include(CppCheck)

# ---- More Options ----

# Profiling
if((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
    option(OGS_PROFILE
           "Enables compiling with flags set for profiling with gprof." OFF)
endif() # GCC AND GPROF_PATH

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

option(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)
option(OGS_INSITU "Builds OGS with insitu visualization capabilities." OFF)
option(OGS_USE_LIS "Use Lis" OFF)
option(OGS_USE_PETSC "Use PETSc routines" OFF)
option(OGS_USE_NETCDF "Add NetCDF support." OFF)

# Eigen
option(OGS_USE_EIGEN "Use Eigen linear solver" ON)
option(OGS_USE_EIGEN_UNSUPPORTED "Use Eigen unsupported modules" ON)
option(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN "" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# We assume that it's save to use vectorization with Eigen < 3.3 (strictly
# smaller than 3.3.!). At least we only observed vectorization issues with Eigen
# 3.3.x. If you want to use Eigen vectorization, make sure that you run all the
# ctests several times, e.g.: $ ctest --repeat-until-fail 50 You might also want
# to take a look at https://github.com/ufz/ogs/issues/1881.
option(EIGEN_DONT_VECTORIZE "Disables explicit vectorization when defined." ON)

set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
    "Default"
    CACHE STRING "Use dynamically allocated shape matrices")
set_property(CACHE OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
             PROPERTY STRINGS "Default" "ON" "OFF")

if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES STREQUAL "Default")
    if(CMAKE_BUILD_TYPE STREQUAL "Release"
       OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL OFF)
    else()
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL ON)
    endif()
else()
    set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL
        ${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES})
endif()
# Eigen End

option(OGS_USE_MKL "Use Intel MKL" OFF)

# Debug
option(OGS_FATAL_ABORT "Abort in OGS_FATAL" OFF)

# Compiler flags
set(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")
option(STL_NO_DEBUG "Disable STL debug in debug build" OFF)

# Print CMake variable values
if(OGS_CMAKE_DEBUG)
    include(ListAllCMakeVariableValues)
    list_all_cmake_variable_values()
endif()

# Code coverage
option(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)

# Options controlling which FEM elements will be compiled
set(OGS_MAX_ELEMENT_DIM
    3
    CACHE STRING "Maximum dimension of FEM elements to be built.")
set(OGS_MAX_ELEMENT_ORDER
    2
    CACHE STRING "Maximum order of FEM elements to be built.")
option(OGS_ENABLE_ELEMENT_SIMPLEX
       "Build FEM elements for simplices (triangles, tetrahedra)." ON)
option(OGS_ENABLE_ELEMENT_CUBOID
       "Build FEM elements for cuboids (quads, hexahedra)." ON)
option(OGS_ENABLE_ELEMENT_PRISM "Build FEM elements for prisms." ON)
option(OGS_ENABLE_ELEMENT_PYRAMID "Build FEM elements for pyramids." ON)

option(OGS_CHECK_HEADER_COMPILATION "Check header for standalone compilation."
       OFF)

option(OGS_USE_PYTHON "Interface with Python" OFF)
if(OGS_USE_PYTHON AND NOT Python3_FOUND)
    message(FATAL_ERROR "Python not found but required for OGS_USE_PYTHON=ON!")
endif()

option(OGS_USE_MFRONT
       "Enable solid material models by MFront (https://tfel.sourceforge.net)"
       OFF)
if(OGS_USE_MFRONT)
    add_definitions(-DOGS_USE_MFRONT)
endif()
# ---- Definitions ----
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
endif()

if(OGS_USE_MPI)
    add_definitions(-DUSE_MPI)
endif()

# Eigen
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()

if(EIGEN_NO_DEBUG)
    add_definitions(-DEIGEN_NO_DEBUG)
endif()

if(EIGEN_DONT_VECTORIZE)
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
# End Eigen

if(OGS_FATAL_ABORT)
    add_definitions(-DOGS_FATAL_ABORT)
endif()

# Packaging
include(scripts/cmake/packaging/Pack.cmake)

# ---- Subdirectories ----

# External projects
if(OGS_INSITU)
    include("${PARAVIEW_USE_FILE}")
    add_definitions(-DUSE_INSITU)
endif()

add_subdirectory(ThirdParty)
include_directories(
    SYSTEM
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/googletest/googletest/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/autocheck/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/json/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/exprtk
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/spdlog/include)

# vtkdiff
if(EXISTS ${PROJECT_SOURCE_DIR}/ThirdParty/vtkdiff/CMakeLists.txt AND BUILD_TESTING)
    add_subdirectory(ThirdParty/vtkdiff)
    install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin COMPONENT ogs_extras)
endif()

include(scripts/cmake/CheckHeaderCompilation.cmake)

add_subdirectory(Applications)
add_subdirectory(BaseLib)
add_subdirectory(GeoLib)
add_subdirectory(InfoLib)
add_subdirectory(MathLib)
add_subdirectory(MeshLib)
add_subdirectory(MeshGeoToolsLib)
add_subdirectory(NumLib)

if(OGS_BUILD_PROCESS_ComponentTransport
   OR OGS_BUILD_PROCESS_RichardsComponentTransport)
    add_subdirectory(ChemistryLib)
endif()

if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR BUILD_TESTING)
    add_subdirectory(ParameterLib)
    add_subdirectory(MaterialLib)
    add_subdirectory(ProcessLib)
endif()
if(BUILD_TESTING AND NOT IS_SUBPROJECT)
    add_subdirectory(Tests)
endif()

include(UnityBuildSettings)

file(WRITE ${PROJECT_BINARY_DIR}/disabled-tests.log "${DISABLED_TESTS_LOG}")
unset(DISABLED_TESTS_LOG CACHE) # Don't write to CMakeCache.txt

check_header_compilation()

include(scripts/cmake/MarkVariablesAdvanced.cmake)
unset(PRE_INSTALL_RUN CACHE)

include(FeatureSummary)
feature_summary(WHAT ALL)
