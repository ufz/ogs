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
if(OGS_USE_MFRONT)
    set(CMAKE_REQUIRE_FIND_PACKAGE_Threads TRUE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif()
find_package(Threads)

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

if(OGS_USE_NETCDF)
    set(NETCDF_ROOT ${CONAN_NETCDF-C_ROOT})
    set(NETCDF_CXX_ROOT ${CONAN_NETCDF-CXX_ROOT})
    find_package(NetCDF REQUIRED)
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

# blas / lapack
if(OGS_USE_MKL)
    if("${MKL_USE_interface}" STREQUAL "lp64")
        set(BLA_VENDOR Intel10_64lp)
    elseif("${MKL_USE_interface}" STREQUAL "ilp64")
        set(BLA_VENDOR Intel10_64ilp)
    endif()
endif()
find_package(BLAS)
find_package(LAPACK)

if(OGS_USE_MKL)
    find_package(MKL REQUIRED)
    find_file(MKL_SETVARS setvars.sh PATHS ${MKL_ROOT_DIR} ${MKL_ROOT_DIR}/..
                                           ${MKL_ROOT_DIR}/../..
              NO_DEFAULT_PATH
    )
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
                    "NOTE: Please add the MKL redist directory to your PATH environment variable!\nE.g. with: set PATH=%PATH%;${MKL_ROOT_DIR}/redist/intel64"
            )
        else()
            message(
                STATUS
                    "NOTE: Please set LD_LIBRARY_PATH with `export LD_LIBRARY_PATH=${MKL_LIBRARY_DIR}`!\n"
            )
        endif()
    endif()
endfunction()
