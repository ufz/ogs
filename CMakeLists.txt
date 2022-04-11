# ---- OGS-6 Project ----
cmake_minimum_required(VERSION 3.22)

project(OGS-6)

include(scripts/cmake/CPM.cmake)
include(scripts/cmake/CMakeSetup.cmake)
include(ParseCMakeArgs)
if(DEFINED ENV{CI})
    message(STATUS "CMake arguments: ${CMAKE_ARGS}")
endif()
include(CTest)

# ---- Preliminary Options ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "Create shared libraries?" OFF)
option(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
set(CMAKE_LIBRARY_SEARCH_PATH ""
    CACHE PATH
          "Additional library installation path, e.g. /opt/local or C:/libs"
)
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, \
    defaults to native (*nix) / blend (MSVC)."
)
option(OGS_USE_CONAN "Should Conan package manager be used?" ON)
set(OGS_CONAN_BUILD "missing" CACHE STRING "Possible values: all, missing, \
    never or list of libs to build"
)
option(OGS_DISABLE_COMPILER_CACHE "Disables compiler cache." OFF)
option(OGS_USE_UNITY_BUILDS "Enables Unity builds for faster compilation." ON)
option(OGS_USE_PYTHON "Interface with Python" ON)
option(OGS_USE_PIP "Enables automatic Python virtual environment handling." OFF)
cmake_dependent_option(
    OGS_BUILD_SWMM "Should the SWMM interface be built?" ON "WIN32" OFF
)

if(OGS_USE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()
set(OGS_PETSC_CONFIG_OPTIONS "" CACHE STRING
                                      "Additional PETSc configuration options."
)
option(
    OGS_INSTALL_EXTERNAL_DEPENDENCIES
    "Should the external dependencies (TFEL, PETSc, LIS) be installed? Set CMAKE_INSTALL_PREFIX as well!"
    OFF
)
option(OGS_USE_CVODE "Use the Sundials CVODE module?" OFF)
option(OGS_BUILD_UTILS "Should the utilities programs be built?" ON)
option(OGS_BUILD_TESTING "Should the tests be built?" ON)

if(MSVC)
    set(CMD_COMMAND "cmd;/c")
endif()

option(OGS_USE_MKL "Use Intel MKL" OFF)
if(OGS_USE_MKL)
    option(MKL_USE_parallel "Use MKL parallel" True)
    option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
    set(MKL_USE_interface
        "lp64"
        CACHE
            STRING
            "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall"
    )
endif()

# ---- CMake includes ----
include(Versions)
include(GitSetup)
include(PythonSetup)
include(ProcessesSetup)
include(ProjectSetup)
include(Functions)
include(ConanSetup)
include(CompilerSetup)
include(JobPools)
include(Find)
include(packaging/PackagingDependencies)
include(DependenciesExternalProject)
include(Dependencies)
include(DocumentationSetup)
include(test/Test)
if(OGS_COVERAGE AND NOT _IS_SUBPROJECT)
    include(Coverage)
endif()
include(CppCheck)

# ---- More Options ----

# Profiling
if((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
    option(OGS_PROFILE
           "Enables compiling with flags set for profiling with gprof." OFF
    )
endif() # GCC AND GPROF_PATH

option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
option(OGS_INSITU "Builds OGS with insitu visualization capabilities." OFF)
option(OGS_USE_LIS "Use Lis" OFF)
option(OGS_USE_PETSC "Use PETSc routines" OFF)
if(OGS_USE_PETSC AND MSVC)
    message(
        FATAL_ERROR
            "OGS_USE_PETSC=ON is not supported on Windows Visual Studio! Use Linux or macOS."
    )
endif()
option(OGS_USE_NETCDF "Add NetCDF support." OFF)

# Eigen
option(OGS_USE_EIGEN_UNSUPPORTED "Use Eigen unsupported modules" ON)
option(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN "" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# We assume that it's save to use vectorization with Eigen < 3.3 (strictly
# smaller than 3.3.!). At least we only observed vectorization issues with Eigen
# 3.3.x. If you want to use Eigen vectorization, make sure that you run all the
# ctests several times, e.g.: $ ctest --repeat-until-fail 50 You might also want
# to take a look at https://github.com/ufz/ogs/issues/1881.
option(EIGEN_DONT_VECTORIZE "Disables explicit vectorization when defined." ON)

set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Default"
    CACHE STRING "Use dynamically allocated shape matrices"
)
set_property(
    CACHE OGS_EIGEN_DYNAMIC_SHAPE_MATRICES PROPERTY STRINGS "Default" "ON"
                                                    "OFF"
)

if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES STREQUAL "Default")
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL
                                              "RelWithDebInfo"
    )
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL OFF)
    else()
        set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL ON)
    endif()
else()
    set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL
        ${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES}
    )
endif()

cmake_dependent_option(
    OGS_EIGEN_USE_MKL "Defines EIGEN_USE_MKL_ALL." OFF "OGS_USE_MKL" OFF
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

include(packaging/Pack)

# ---- Subdirectories ----
include_directories(${PROJECT_SOURCE_DIR})
# xdmfdiff
if(OGS_BUILD_TESTING AND TARGET OgsXdmf)
    add_subdirectory(Tests/xdmfdiff)
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

if(_build_chemistry_lib)
    add_subdirectory(ChemistryLib)
endif()

if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR OGS_BUILD_TESTING)
    add_subdirectory(ParameterLib)
    add_subdirectory(MaterialLib)
    add_subdirectory(ProcessLib)
endif()
if(OGS_BUILD_TESTING AND NOT _IS_SUBPROJECT)
    add_subdirectory(Tests)
endif()

include(UnityBuildSettings)

if(OGS_USE_PIP)
    # Caches a hash of requested Python packages when they were successfully
    # installed. On subsequent runs compare new hash to cached. If equal do
    # nothing.
    get_property(
        _addtest_python_packages GLOBAL PROPERTY AddTest_PYTHON_PACKAGES
    )
    list(APPEND OGS_PYTHON_PACKAGES ${_addtest_python_packages})
    list(REMOVE_DUPLICATES OGS_PYTHON_PACKAGES)
    list(SORT OGS_PYTHON_PACKAGES)
    string(SHA1 _ogs_python_packages_sha1 "${OGS_PYTHON_PACKAGES}")
    list(LENGTH OGS_PYTHON_PACKAGES OGS_PYTHON_PACKAGES_LENGTH)
    if(NOT ${_ogs_python_packages_sha1} STREQUAL "${_OGS_PYTHON_PACKAGES_SHA1}"
       AND ${OGS_PYTHON_PACKAGES_LENGTH} GREATER 0
    )
        string(REPLACE ";" "\n" REQUIREMENTS_CONTENT "${OGS_PYTHON_PACKAGES}")
        file(WRITE ${PROJECT_BINARY_DIR}/requirements.txt
             ${REQUIREMENTS_CONTENT}
        )
        execute_process(
            COMMAND ${LOCAL_VIRTUALENV_BIN_DIR}/pip install -r requirements.txt
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            RESULT_VARIABLE _return_code
            OUTPUT_VARIABLE _out
            ERROR_VARIABLE _err
        )
        if(${_return_code} EQUAL 0)
            set(_OGS_PYTHON_PACKAGES_SHA1 "${_ogs_python_packages_sha1}"
                CACHE INTERNAL ""
            )
        else()
            message(
                FATAL_ERROR
                    "Installation of Python packages via pip failed!\n"
                    "To disable pip set OGS_USE_PIP=OFF.\n\n${_out}\n${_err}"
            )
        endif()
    endif()
endif()

file(WRITE ${PROJECT_BINARY_DIR}/disabled-tests.log "${DISABLED_TESTS_LOG}")
unset(DISABLED_TESTS_LOG CACHE) # Don't write to CMakeCache.txt

check_header_compilation()

include(scripts/cmake/MarkVariablesAdvanced.cmake)
unset(PRE_INSTALL_RUN CACHE)

include(Features)

printEnabledProcesses()

printMKLUsage()
