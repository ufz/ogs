if(BUILD_TESTING)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        GIT_TAG 389cb68b87193358358ae87cc56d257fd0d80189
        OPTIONS
            "INSTALL_GTEST OFF"
            "gtest_force_shared_crt ON"
    )

    CPMAddPackage(
        NAME autocheck
        GITHUB_REPOSITORY ufz/autocheck
        GIT_TAG e388ecbb31c49fc2724c8d0436da313b6edca7fd
        DOWNLOAD_ONLY YES
    )
    if(autocheck_ADDED)
        add_library(autocheck INTERFACE IMPORTED)
        target_include_directories(autocheck INTERFACE ${autocheck_SOURCE_DIR}/include)
    endif()

    CPMAddPackage(
        NAME vtkdiff
        GITHUB_REPOSITORY ufz/vtkdiff
        GIT_TAG 49403cee266bb8e80405a02d677dbb5f71afc61a
        OPTIONS
            "VTK_LIBRARIES vtkIOXML"
    )
    if(vtkdiff_ADDED)
        install(PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin)
    endif()
endif()

CPMAddPackage(
    NAME exprtk
    GITHUB_REPOSITORY ArashPartow/exprtk
    GIT_TAG c7c219480d9678eec7383a4a99030683c4a84d91
    DOWNLOAD_ONLY YES
)
if(exprtk_ADDED)
    add_library(exprtk INTERFACE IMPORTED)
    target_include_directories(exprtk INTERFACE ${exprtk_SOURCE_DIR})
endif()

CPMAddPackage(
    NAME spdlog
    GITHUB_REPOSITORY gabime/spdlog
    VERSION 1.8.2
)

CPMAddPackage(
    NAME tclap
    GITHUB_REPOSITORY ufz/tclap
    GIT_TAG 03abc3a3327214137c6ffd5b9a6efe23f0927cc2
    DOWNLOAD_ONLY YES
)
if(tclap_ADDED)
    add_library(tclap INTERFACE IMPORTED)
    target_include_directories(tclap INTERFACE ${tclap_SOURCE_DIR}/include)
endif()

CPMAddPackage(
    NAME tetgen
    GITHUB_REPOSITORY ufz/tetgen
    GIT_TAG 603ba181ebfaed38eec88532720e282606009b73
)
if(tetgen_ADDED)
    install(PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin)
endif()

######################
### Find tools     ###
######################

string(REPLACE ".windows.1" "" GIT_VERSION_STRING ${GIT_VERSION_STRING})
if(${GIT_VERSION_STRING} VERSION_LESS ${ogs.minimum_version.git})
    message(FATAL_ERROR "Git version ${ogs.minimum_version.git} is required. \
        Found version ${GIT_VERSION_STRING}.")
endif()

find_package(Doxygen OPTIONAL_COMPONENTS dot)

# Find gnu profiler gprof
find_program(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)

find_program(CPPCHECK_TOOL_PATH cppcheck)

# Find bash itself ...
find_program(BASH_TOOL_PATH bash
    HINTS ${GITHUB_BIN_DIR} DOC "The bash executable")

# Dumpbin is a windows dependency analaysis tool required for packaging.
# Variable has to be named gp_cmd to override the outdated find routines
# of the GetPrerequisites CMake-module.
if(WIN32)
    include(MSVCPaths)
    find_program(gp_cmd dumpbin DOC "Windows dependency analysis tool"
        PATHS ${MSVC_INSTALL_PATHS} PATH_SUFFIXES VC/bin)
    if(gp_cmd)
        get_filename_component(dir ${gp_cmd} PATH)
        set(ENV{PATH} "${dir}/../../../Common7/IDE;$ENV{PATH}")
    endif()
endif()

find_program(CURL_TOOL_PATH curl DOC "The curl-tool")

find_program(S3CMD_TOOL_PATH s3cmd DOC "S3cmd tool for uploading to Amazon S3")

if(WIN32)
    find_program(CLCACHE_TOOL_PATH clcache)
else()
    find_program(CCACHE_TOOL_PATH ccache)
endif()

# Tools for web
find_program(VTKJS_CONVERTER vtkDataConverter
    PATHS ${PROJECT_SOURCE_DIR}/web/node_modules/.bin)
find_program(HUGO hugo)
find_program(NPM npm)
find_program(YARN yarn)
find_program(PIP pip)
find_program(PANDOC_CITEPROC pandoc-citeproc)

find_program(MODULE_CMD lmod PATHS /software/lmod/lmod/libexec)

find_program(SNAKEMAKE snakemake HINTS ${LOCAL_VIRTUALENV_BIN_DIRS})
find_program(PARSL parsl-visualize HINTS ${LOCAL_VIRTUALENV_BIN_DIRS})

find_program(GMSH gmsh)

######################
### Find libraries ###
######################
find_package(Boost ${ogs.minimum_version.boost} REQUIRED)

set(VTK_COMPONENTS vtkIOXML vtkIOLegacy)
if(OGS_BUILD_GUI)
    set(VTK_COMPONENTS ${VTK_COMPONENTS}
        vtkIOExport vtkImagingCore
        vtkInteractionStyle vtkInteractionWidgets
        vtkGUISupportQt vtkRenderingOpenGL2 vtkRenderingContextOpenGL2
        vtkFiltersTexture vtkRenderingAnnotation vtkRenderingCore vtkFiltersParallel
    )
endif()
if(OGS_USE_MPI)
    set(VTK_COMPONENTS ${VTK_COMPONENTS} vtkIOParallelXML vtkParallelMPI)
endif()
if(OGS_INSITU)
    find_package(ParaView REQUIRED)
else()
    find_package(VTK ${ogs.minimum_version.vtk} REQUIRED COMPONENTS ${VTK_COMPONENTS})
    include(${VTK_USE_FILE})
endif()

if(OGS_USE_CONAN)
    set(EIGEN3_INCLUDE_DIR ${CONAN_INCLUDE_DIRS_EIGEN} CACHE INTERNAL "")
else()
    find_package(Eigen3 ${ogs.minimum_version.eigen} REQUIRED)
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

if(OGS_USE_MFRONT)
    ## pthread, is a requirement of mfront ##
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    if(CMAKE_USE_PTHREADS_INIT)
        set(HAVE_PTHREADS TRUE)
        add_definitions(-DHAVE_PTHREADS)
    endif()
    if(OGS_USE_CONAN)
        set(TFELHOME ${CONAN_TFEL_ROOT} CACHE INTERNAL "")
    endif()
endif()

# Do not search for libs if this option is set
if(OGS_NO_EXTERNAL_LIBS)
    return()
endif() # OGS_NO_EXTERNAL_LIBS

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## Qt5 library ##
if(OGS_BUILD_GUI)
    set(QT_MODULES Gui Widgets Xml XmlPatterns)
    if(UNIX AND NOT APPLE)
        list(APPEND QT_MODULES X11Extras)
    endif()
    find_package(Qt5 ${ogs.minimum_version.qt} REQUIRED ${QT_MODULES})
    cmake_policy(SET CMP0020 NEW)
endif()

if(OGS_USE_NETCDF)
    set(NETCDF_ROOT ${CONAN_NETCDF-C_ROOT})
    set(NETCDF_CXX_ROOT ${CONAN_NETCDF-CXX_ROOT})
    find_package(NetCDF REQUIRED)
    add_compile_options(-DOGS_USE_NETCDF)
endif()

# lapack
find_package(LAPACK QUIET)

## geotiff ##
find_package(GEOTIFF)

## lis ##
if(OGS_USE_LIS)
    find_package( LIS REQUIRED )
endif()

if(OGS_USE_MKL)
    find_package( MKL REQUIRED )
endif()

if(OGS_USE_PETSC)
    message(STATUS "Configuring for PETSc")

    option(FORCE_PETSC_EXECUTABLE_RUNS
        "Force CMake to accept a given PETSc configuration" ON)

    # Force CMake to accept a given PETSc configuration in case the failure of
    # MPI tests. This may cause the compilation broken.
    if(FORCE_PETSC_EXECUTABLE_RUNS)
        set(PETSC_EXECUTABLE_RUNS YES)
    endif()

    find_package(PETSc ${ogs.minimum_version.petsc} REQUIRED)

    include_directories(SYSTEM ${PETSC_INCLUDES})
endif()

## Check MPI package
if(OGS_USE_MPI)
    find_package(MPI REQUIRED)
endif()

find_package(Shapelib)
if(Shapelib_FOUND)
    include_directories(SYSTEM ${Shapelib_INCLUDE_DIRS})
elseif(OGS_BUILD_GUI)
    message(FATAL_ERROR "Shapelib not found but it is required for OGS_BUILD_GUI!")
endif()

## Sundials cvode ode-solver library
if(OGS_USE_CVODE)
    find_package(CVODE REQUIRED)
endif()

find_package(Filesystem REQUIRED COMPONENTS Final Experimental)
configure_file(${PROJECT_SOURCE_DIR}/BaseLib/filesystem.h.in
               ${PROJECT_BINARY_DIR}/BaseLib/filesystem.h)
