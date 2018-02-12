######################
### Find tools     ###
######################

# Find doxygen
if(WIN32)
    find_program(DOXYGEN_DOT_EXECUTABLE NAMES dot
        PATHS "$ENV{ProgramFiles}/Graphviz*/bin")
    find_package(Doxygen QUIET)
    if(DOXYGEN_DOT_PATH)
        file(TO_NATIVE_PATH ${DOXYGEN_DOT_PATH} DOXYGEN_DOT_PATH)
        set(DOXYGEN_DOT_PATH "\"${DOXYGEN_DOT_PATH}\"")
    endif()
else()
    find_package(Doxygen QUIET)
endif()

# Find gnu profiler gprof
find_program(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)

find_package(cppcheck QUIET)

find_package(PythonInterp QUIET)

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

find_program(CCACHE_TOOL_PATH ccache)

# Tools for web
find_program(VTKJS_CONVERTER vtkDataConverter
    PATHS ${PROJECT_SOURCE_DIR}/web/node_modules/.bin)
find_program(HUGO hugo)
find_program(NPM npm)
find_program(YARN yarn)
find_program(PIP pip)
find_package(PythonInterp)
find_program(PANDOC_CITEPROC pandoc-citeproc)

find_program(MODULE_CMD modulecmd
    PATHS /usr/local/modules/3.2.10-1/Modules/3.2.10/bin)

######################
### Find libraries ###
######################
find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include(VtkModules)
find_package(VTK COMPONENTS ${VTK_MODULES} REQUIRED)
include(${VTK_USE_FILE})
if(NOT OGS_USE_CONAN)
    foreach(DIR ${VTK_INCLUDE_DIRS})
        if("${DIR}" MATCHES ".*vtknetcdf.*")
            include_directories(SYSTEM ${DIR}/../cxx ${DIR}/include)
        elseif("${DIR}" MATCHES ".*vtk.*")
            include_directories(SYSTEM ${DIR}/vtknetcdf/include)
        endif()
    endforeach()
    include_directories(SYSTEM ${VTK_DIR}/../ThirdParty/netcdf/vtknetcdf/cxx)
endif()

find_package(Eigen3 3.2.9 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

## pthread, is a requirement of logog ##
if(CMAKE_CROSSCOMPILING)
    set(THREADS_PTHREAD_ARG 0 CACHE STRING "Result from TRY_RUN" FORCE)
endif()
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if(CMAKE_USE_PTHREADS_INIT)
    set(HAVE_PTHREADS TRUE)
    add_definitions(-DHAVE_PTHREADS)
endif()

# Do not search for libs if this option is set
if(OGS_NO_EXTERNAL_LIBS)
    return()
endif() # OGS_NO_EXTERNAL_LIBS

find_package(OpenMP QUIET)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP enabled.")
endif()

find_package(Metis QUIET)

## Qt5 library ##
if(OGS_BUILD_GUI)
    find_package(Qt5 5.2 REQUIRED Gui Widgets Xml XmlPatterns)
    cmake_policy(SET CMP0020 NEW)
    set(CMAKE_AUTOMOC TRUE)
    set(CMAKE_AUTOUIC TRUE)
    set(CMAKE_AUTORCC TRUE)
    if(CMAKE_CROSSCOMPILING)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(QT_XML_DEPS REQUIRED Xml)
        list(REMOVE_ITEM QT_XML_DEPS_LIBRARIES Xml Core)
        pkg_check_modules(QT_GUI_DEPS REQUIRED Gui)
        list(REMOVE_ITEM QT_GUI_DEPS_LIBRARIES Gui Core)
        pkg_check_modules(QT_NETWORK_DEPS REQUIRED Network)
        list(REMOVE_ITEM QT_NETWORK_DEPS_LIBRARIES Network Core)
    endif()
endif()

# lapack
find_package(LAPACK QUIET)

## geotiff ##
find_package(LibGeoTiff)
if(GEOTIFF_FOUND)
    add_definitions(-DGEOTIFF_FOUND)
endif() # GEOTIFF_FOUND

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

    find_package(PETSc REQUIRED)

    include_directories(SYSTEM ${PETSC_INCLUDES})

    add_definitions(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)

endif()

find_package(OpenSSL)

## Check MPI package
if(OGS_USE_MPI)
    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
endif()

find_package(Shapelib)
if(Shapelib_FOUND)
    include_directories(SYSTEM ${Shapelib_INCLUDE_DIRS})
elseif(OGS_BUILD_GUI)
    message(FATAL_ERROR "Shapelib not found but it is required for OGS_BUILD_GUI!")
endif()

## Sundials cvode ode-solver library
find_package(CVODE)
if(CVODE_FOUND)
    add_definitions(-DCVODE_FOUND)
endif() # CVODE_FOUND
