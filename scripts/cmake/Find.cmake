# ##############################################################################
# Find tools
# ##############################################################################

find_package(Doxygen 1.9.2 OPTIONAL_COMPONENTS dot)
if(TARGET Doxygen::dot)
    # Create dependency graph in build dir with:
    # ~~~
    # cmake . --graphviz=cmake-dependencies.dot && \
    #   dot -Tpng cmake-dependencies.dot -o cmake-dependencies.png
    # ~~~
    file(WRITE ${PROJECT_BINARY_DIR}/CMakeGraphVizOptions.cmake
         "set(GRAPHVIZ_IGNORE_TARGETS testrunner \"vtk.*\")\n"
    )
endif()

# Find gnu profiler gprof
find_program(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)

find_program(CPPCHECK_TOOL_PATH cppcheck)

# Find bash itself ...
find_program(BASH_TOOL_PATH bash DOC "The bash executable")

find_program(CCACHE_TOOL_PATH ccache)

find_program(MODULE_CMD lmod PATHS /software/lmod/lmod/libexec)

find_program(SNAKEMAKE snakemake)
find_program(TEE_TOOL_PATH tee)
if(OGS_BUILD_TESTING AND SNAKEMAKE AND NOT TEE_TOOL_PATH)
    message(WARNING "tee tool was not found. Snakemake tests are disabled!")
endif()

find_program(GMSH gmsh)

find_program(XMLSTARLET_TOOL_PATH xmlstarlet)

if(OGS_INCLUDE_WHAT_YOU_USE)
    find_program(INCLUDE_WHAT_YOU_USE_TOOL_PATH include-what-you-use REQUIRED)
endif()

# ##############################################################################
# Find libraries
# ##############################################################################

# Qt5 library ##
if(OGS_BUILD_GUI)
    set(QT_MODULES Gui Widgets Xml XmlPatterns)
    if(UNIX AND NOT APPLE)
        list(APPEND QT_MODULES X11Extras)
    endif()
    find_package(Qt5 ${ogs.minimum_version.qt} REQUIRED ${QT_MODULES})
    cmake_policy(SET CMP0020 NEW)
    list(APPEND CMAKE_INSTALL_RPATH ${Qt5_DIR}/../..)
    list(APPEND CMAKE_BUILD_RPATH ${Qt5_DIR}/../..)
endif()

# geotiff ##
find_package(GEOTIFF)

cmake_dependent_option(
    OGS_USE_PETSC_MKL_EIGEN_OPENMP
    "When PETSc and MKL is used, shall OpenMP be used for Eigen (or Intels iomp if false (default))?"
    OFF
    "OGS_USE_PETSC;OGS_USE_MKL"
    OFF
)
if(NOT (OGS_USE_PETSC AND OGS_USE_MKL) OR OGS_USE_PETSC_MKL_EIGEN_OPENMP)
    # this pulls in libgomp dependency, when MKL is enabled libiomp5 is used.
    find_package(OpenMP COMPONENTS C CXX)
endif()

# blas / lapack / MKL
if(OGS_USE_MKL)
    if(APPLE)
        set(_mac_ld_prefix "DY")
    endif()
    if(NOT GUIX_BUILD AND (NOT DEFINED ENV{MKLROOT} OR
                          (NOT "$ENV{${_mac_ld_prefix}LD_LIBRARY_PATH}"
                                    MATCHES "intel" AND NOT WIN32))
    )
        message(
            FATAL_ERROR
                "OGS_USE_MKL was used but it seems that you did not source the MKL environment. "
                "Typically you can run `source /opt/intel/oneapi/setvars.sh` before running CMake."
        )
    endif()
    set(MKL_INTERFACE
        "lp64"
        CACHE
            STRING
            "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall"
    )
    if("${MKL_INTERFACE}" STREQUAL "lp64")
        set(BLA_VENDOR Intel10_64lp)
    elseif("${MKL_INTERFACE}" STREQUAL "ilp64")
        set(BLA_VENDOR Intel10_64ilp)
    endif()
    if(NOT WIN32 AND NOT APPLE AND NOT GUIX_BUILD)
        set(CMAKE_REQUIRE_FIND_PACKAGE_BLAS TRUE)
        set(CMAKE_REQUIRE_FIND_PACKAGE_LAPACK TRUE)
    endif()
endif()
find_package(BLAS)
find_package(LAPACK)

if(OGS_USE_MKL)
    if(GUIX_BUILD)
        find_package(PkgConfig REQUIRED)
        if(DEFINED ENV{GUIX_ENVIRONMENT})
            set(MKLROOT $ENV{GUIX_ENVIRONMENT})
        endif()
        set(PKG_CONFIG_ARGN "--define-variable=MKLROOT=${MKLROOT}")
        set(ENV{PKG_CONFIG_PATH} ${MKLROOT}/bin/pkgconfig)
        # TODO: Using -seq instead of -iomp; library iomp5 is missing
        # See https://gitlab.opengeosys.org/ogs/inf/guix-ogs/-/issues/1
        pkg_search_module(MKL REQUIRED IMPORTED_TARGET mkl-dynamic-${MKL_INTERFACE}-seq)
        add_library(MKL::MKL ALIAS PkgConfig::MKL)
    else()
        find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT} ${MKLROOT})
        find_file(MKL_SETVARS setvars.sh PATHS ${MKL_ROOT}/../.. NO_DEFAULT_PATH)
    endif()
endif()

# Check MPI package
if(OGS_USE_MPI)
    find_package(MPI REQUIRED)
    if(NOT "${MPI_C_COMPILER}" STREQUAL "${CMAKE_C_COMPILER}"
       OR NOT "${MPI_CXX_COMPILER}" STREQUAL "${CMAKE_CXX_COMPILER}"
    )
        message(
            FATAL_ERROR
                "The selected compilers\n"
                "  - ${CMAKE_C_COMPILER} and \n"
                "  - ${CMAKE_CXX_COMPILER}\n"
                "are not MPI-enabled compiler!\n"
                "Set compiler on a clean build directory with e.g.:\n"
                "export CC=`which mpicc` CXX=`which mpic++`\n"
                "cmake ../ogs -DOGS_USE_PETSC=ON"
        )
    endif()
endif()

# Prints instructions for setting MKL runtime environment.
function(printMKLUsage)
    if(NOT OGS_USE_MKL)
        return()
    endif()
    if(MKL_SETVARS)
        message(
            STATUS
                "NOTE: Please run `source ${MKL_SETVARS}` to set LD_LIBRARY_PATH for MKL!\n"
        )
    else()
        if(WIN32)
            message(
                STATUS
                    "NOTE: In addition to loading the oneAPI setvars.bat file please also add this to your PATH C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\windows\\redist\\intel64_win\\compiler"
            )
        else()
            message(
                STATUS
                    "NOTE: Please set LD_LIBRARY_PATH with `export LD_LIBRARY_PATH=${MKL_LIBRARY_DIR}`!\n"
            )
        endif()
    endif()
endfunction()
