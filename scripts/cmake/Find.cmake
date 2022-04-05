# ##############################################################################
# Find tools
# ##############################################################################

find_package(Doxygen OPTIONAL_COMPONENTS dot)
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
find_program(
    BASH_TOOL_PATH bash HINTS ${GITHUB_BIN_DIR} DOC "The bash executable"
)

find_program(CCACHE_TOOL_PATH ccache)

# Tools for web
find_program(
    VTKJS_CONVERTER vtkDataConverter
    PATHS ${PROJECT_SOURCE_DIR}/web/node_modules/.bin
)

find_program(MODULE_CMD lmod PATHS /software/lmod/lmod/libexec)

find_program(SNAKEMAKE snakemake)
find_program(TEE_TOOL_PATH tee)
if(OGS_BUILD_TESTING AND SNAKEMAKE AND NOT TEE_TOOL_PATH)
    message(WARNING "tee tool was not found. Snakemake tests are disabled!")
endif()

find_program(GMSH gmsh)

# ##############################################################################
# Find libraries
# ##############################################################################
if(OGS_USE_MFRONT)
    # pthread, is a requirement of mfront ##
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
        set(HAVE_PTHREADS TRUE)
        add_definitions(-DHAVE_PTHREADS)
    endif()
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}"
    )
endif()

# Qt5 library ##
if(OGS_BUILD_GUI)
    set(QT_MODULES Gui Widgets Xml XmlPatterns)
    if(UNIX AND NOT APPLE)
        list(APPEND QT_MODULES X11Extras)
    endif()
    find_package(Qt5 ${ogs.minimum_version.qt} REQUIRED ${QT_MODULES})
    cmake_policy(SET CMP0020 NEW)
    list(APPEND CMAKE_INSTALL_RPATH ${Qt5_DIR}/../..)
endif()

if(OGS_USE_NETCDF)
    set(NETCDF_ROOT ${CONAN_NETCDF-C_ROOT})
    set(NETCDF_CXX_ROOT ${CONAN_NETCDF-CXX_ROOT})
    find_package(NetCDF REQUIRED)
    add_compile_options(-DOGS_USE_NETCDF)
endif()

# lapack
find_package(LAPACK QUIET)

# geotiff ##
find_package(GEOTIFF)

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
