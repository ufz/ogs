# Public OGS configuration options.

include(CMakeDependentOption)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(OGS_BUILD_WHEEL "Build Python wheels of OGS" OFF)
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
            "OGS_USE_PETSC=ON is not supported on Windows Visual Studio! "
            "Use Linux or macOS."
    )
endif()
if(OGS_USE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()
set(OGS_PETSC_CONFIG_OPTIONS "" CACHE STRING
                                      "Additional PETSc configuration options."
)
option(OGS_BUILD_UTILS "Should the utilities programs be built?" ON)

set(_ogs_build_testing_default ON)
if(NOT DEFINED OGS_BUILD_TESTING AND DEFINED BUILD_TESTING)
    set(_ogs_build_testing_default "${BUILD_TESTING}")
endif()
option(OGS_BUILD_TESTING "Should the OGS tests be built?"
       ${_ogs_build_testing_default}
)
unset(_ogs_build_testing_default)

option(OGS_USE_MKL "Use Intel MKL" OFF)

# Eigen
option(OGS_USE_EIGEN_UNSUPPORTED "Use Eigen unsupported modules" ON)
option(OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN "" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)
option(EIGEN_DONT_VECTORIZE "Disables explicit vectorization when defined." ON)
set(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Default"
    CACHE STRING "Use dynamically allocated shape matrices"
)
set_property(
    CACHE OGS_EIGEN_DYNAMIC_SHAPE_MATRICES PROPERTY STRINGS "Default" "ON"
                                                    "OFF"
)

# Code coverage
option(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)
option(OGS_COVERAGE_PYTHON "Enables code coverage measurements in Python code."
       OFF
)

# Profiling Find gnu profiler gprof. CMAKE_CXX_COMPILER_ID is available after
# project(), whereas the COMPILER_IS_GCC helper is established later.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    find_program(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)
    if(GPROF_PATH)
        option(OGS_PROFILE
               "Enables compiling with flags set for profiling with gprof." OFF
        )
    endif()
endif()

option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
option(OGS_USE_INSITU "Builds OGS with insitu visualization capabilities." OFF)
option(OGS_USE_LIS "Use Lis" OFF)
option(OGS_USE_NETCDF "Add NetCDF support." OFF)

# Options controlling which FEM elements will be compiled.
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

option(OGS_USE_MFRONT
       "Enable solid material models by MFront (https://tfel.sourceforge.net)"
       OFF
)
option(OGS_INCLUDE_WHAT_YOU_USE "Enable include-what-you-use checks." OFF)
