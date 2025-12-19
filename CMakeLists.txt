# ---- OGS-6 Project ----
cmake_minimum_required(VERSION 3.26)

if(CMAKE_VERSION GREATER_EQUAL 3.29)
    cmake_policy(SET CMP0156 NEW)
endif()
if(CMAKE_VERSION GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()

if(APPLE AND NOT (DEFINED CMAKE_OSX_DEPLOYMENT_TARGET
                  OR DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
)
    # Setting the deployment target to current macOS major version
    cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
    string(REGEX MATCH "^[0-9]+" _os_release_major "${_os_release}")
    set(CMAKE_OSX_DEPLOYMENT_TARGET ${_os_release_major}.0
        CACHE STRING "Minimum macOS operating system version for deployment"
    )
endif()

project(OGS-6)

option(OGS_BUILD_WHEEL "Build Python wheels of OGS" OFF)

include(CMakeDependentOption)
include(scripts/cmake/UseCompilerCache.cmake)
if(NOT OGS_DISABLE_COMPILER_CACHE)
    useCompilerCache()
endif()
include(scripts/cmake/CPM.cmake)
include(scripts/cmake/CMakeSetup.cmake)
include(ParseCMakeArgs)
if(DEFINED ENV{CI})
    message(STATUS "CMake arguments: ${CMAKE_ARGS}")
endif()
include(CTest)

# ---- Preliminary Options ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "Create shared libraries?" ON)
option(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
set(CMAKE_LIBRARY_SEARCH_PATH ""
    CACHE PATH
          "Additional library installation path, e.g. /opt/local or C:/libs"
)
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, \
    defaults to native (*nix) / blend (MSVC)."
)
option(OGS_DISABLE_COMPILER_CACHE "Disables compiler cache." OFF)
option(OGS_USE_UNITY_BUILDS "Enables Unity builds for faster compilation." ON)
option(OGS_USE_PIP "Enables automatic Python virtual environment handling." OFF)
cmake_dependent_option(
    OGS_BUILD_SWMM "Should the SWMM interface be built?" ON "WIN32" OFF
)
option(OGS_USE_PETSC "Use PETSc routines" OFF)
if(OGS_USE_PETSC AND MSVC)
    message(
        FATAL_ERROR
            "OGS_USE_PETSC=ON is not supported on Windows Visual Studio! Use Linux or macOS."
    )
endif()
if(OGS_USE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()
set(OGS_PETSC_CONFIG_OPTIONS "" CACHE STRING
                                      "Additional PETSc configuration options."
)
option(OGS_BUILD_UTILS "Should the utilities programs be built?" ON)
option(OGS_BUILD_TESTING "Should the tests be built?" ON)

if(MSVC)
    set(CMD_COMMAND "cmd;/c")
endif()

option(OGS_USE_MKL "Use Intel MKL" OFF)

# Eigen
option(OGS_USE_EIGEN_UNSUPPORTED "Use Eigen unsupported modules" ON)
option(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN "" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# ---- CMake includes ----
include(Versions)
include(GitSetup)
include(GuixSetup)
include(PythonSetup)
include(ProcessesSetup)
include(ProjectSetup)
include(Functions)
include(CompilerSetup)
include(JobPools)
include(Find)
include(DependenciesExternalProject)
include(Dependencies)
include(DocumentationSetup)
include(test/Test)
include(Coverage)
include(CppCheck)
include(ClangTidy)

# ---- More Options ----
set(CMAKE_LINK_LIBRARIES_ONLY_TARGETS ON)

# Profiling
if((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
    option(OGS_PROFILE
           "Enables compiling with flags set for profiling with gprof." OFF
    )
endif() # GCC AND GPROF_PATH

option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
option(OGS_USE_INSITU "Builds OGS with insitu visualization capabilities." OFF)
option(OGS_USE_LIS "Use Lis" OFF)
option(OGS_USE_NETCDF "Add NetCDF support." OFF)

# We assume that it's save to use vectorization with Eigen < 3.3 (strictly
# smaller than 3.3.!). At least we only observed vectorization issues with Eigen
# 3.3.x. If you want to use Eigen vectorization, make sure that you run all the
# ctests several times, e.g.: $ ctest --repeat-until-fail 50 You might also want
# to take a look at https://gitlab.opengeosys.org/ogs/ogs/-/issues/1881.
option(EIGEN_DONT_VECTORIZE "Disables explicit vectorization when defined." ON)

# When left as "Default":
#  Build type Release or RelWithDebInfo: OFF, otherwise ON, see MathLib/CMakeLists.txt.
set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Default"
    CACHE STRING "Use dynamically allocated shape matrices"
)
set_property(
    CACHE OGS_EIGEN_DYNAMIC_SHAPE_MATRICES PROPERTY STRINGS "Default" "ON"
                                                    "OFF"
)
# Eigen End

# Debug
option(OGS_FATAL_ABORT "Abort in OGS_FATAL" OFF)

# Compiler flags
set(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")

# Print CMake variable values
if(OGS_CMAKE_DEBUG)
    include(ListAllCMakeVariableValues)
    list_all_cmake_variable_values()
endif()

# Code coverage
option(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)

# Options controlling which FEM elements will be compiled
set(OGS_MAX_ELEMENT_DIM 3
    CACHE STRING "Maximum dimension of FEM elements to be built."
)
set(OGS_MAX_ELEMENT_ORDER 2 CACHE STRING
                                  "Maximum order of FEM elements to be built."
)
option(OGS_ENABLE_ELEMENT_SIMPLEX
       "Build FEM elements for simplices (triangles, tetrahedra)." ON
)
option(OGS_ENABLE_ELEMENT_CUBOID
       "Build FEM elements for cuboids (quads, hexahedra)." ON
)
option(OGS_ENABLE_ELEMENT_PRISM "Build FEM elements for prisms." ON)
option(OGS_ENABLE_ELEMENT_PYRAMID "Build FEM elements for pyramids." ON)
if(NOT OGS_MAX_ELEMENT_DIM MATCHES "^[0-3]$")
    message(
        FATAL_ERROR "OGS_MAX_ELEMENT_DIM must be an integer between 0 and 3."
    )
endif()
if(NOT OGS_MAX_ELEMENT_ORDER MATCHES "^[0-9]$")
    message(FATAL_ERROR "OGS_MAX_ELEMENT_ORDER must be an integer.")
endif()

option(OGS_CHECK_HEADER_COMPILATION "Check header for standalone compilation."
       OFF
)

option(OGS_USE_MFRONT
       "Enable solid material models by MFront (https://tfel.sourceforge.net)"
       OFF
)

option(OGS_INCLUDE_WHAT_YOU_USE "Enable include-what-you-use checks." OFF)

include(packaging/Pack)

# ---- Subdirectories ----
include(scripts/cmake/CheckHeaderCompilation.cmake)

add_subdirectory(Applications)
add_subdirectory(BaseLib)
add_subdirectory(GeoLib)
add_subdirectory(InfoLib)
add_subdirectory(MathLib)
add_subdirectory(MeshLib)
add_subdirectory(MeshToolsLib)
add_subdirectory(MeshGeoToolsLib)
add_subdirectory(NumLib)

if(_build_chemistry_lib)
    add_subdirectory(ChemistryLib)
endif()

if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR OGS_BUILD_TESTING)
    if(TARGET OgsXdmf)
        add_subdirectory(Tests/xdmfdiff)
    endif()
    add_subdirectory(ParameterLib)
    add_subdirectory(MaterialLib)
    add_subdirectory(ProcessLib)
endif()
if(OGS_BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
    add_subdirectory(Tests)
endif()

include(UnityBuildSettings)
include(GenerateEnvrcFile)
check_header_compilation()
include(Features)
printEnabledProcesses()
printMKLUsage()

# Check wheel binaries
list(SORT _OGS_EXECUTABLES)
list(REMOVE_ITEM _OGS_EXECUTABLES testrunner SWMMConverter) # SWMMConverter is win-only

if(OGS_BUILD_WHEEL)
    execute_process(
        COMMAND ${Python_EXECUTABLE} scripts/python/compare_binaries.py "${_OGS_EXECUTABLES}"
        ECHO_OUTPUT_VARIABLE
        ECHO_ERROR_VARIABLE
        RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Binary list comparison failed. Please update binaries_list in provide_ogs_cli_tools_via_wheel.py!")
    endif()
endif()

if(TARGET ogs_test_dlopen)
    # Has to come last in this top-level CMakeLists.txt to add tests for all shared libraries.
    # We pass ${Python_LIBRARIES} as additional lib(s) to load, because many OGS libs
    # (e.g. each process, because of a cyclic dependency between the generic process lib and (Python) BCs)
    # depend on them and they are not explicitly linked against.
    for_each_subdir("${PROJECT_SOURCE_DIR}" add_dlopen_tests ${Python_LIBRARIES})
endif()
