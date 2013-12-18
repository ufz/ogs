############################
### Find OGS directories ###
############################

# Compiled libraries (for Windows)
FIND_PATH(OGS_LIBS_DIR_FOUND geotiff.lib
	PATHS $ENV{OGS_LIBS} ${OGS_LIBS_DIR} ${PROJECT_SOURCE_DIR}/../libs C:/OGS_Libs
	PATH_SUFFIXES libgeotiff)
IF(OGS_LIBS_DIR_FOUND)
	SET(OGS_LIBS_DIR ${OGS_LIBS_DIR_FOUND}/.. CACHE STRING "")
ENDIF()

######################
### Find tools     ###
######################

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen QUIET)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)

FIND_PACKAGE(cppcheck QUIET)

FIND_PACKAGE(PythonInterp QUIET)

FIND_PACKAGE(GitHub)

FIND_PROGRAM(GIT_TOOL_PATH git HINTS ${GITHUB_BIN_DIR} DOC "The git command line interface")
IF(NOT GIT_TOOL_PATH)
	IF(WIN32)
		MESSAGE(FATAL_ERROR "Git not found! Please install GitHub for Windows or Git!")
	ELSE()
		MESSAGE(FATAL_ERROR "Git not found but is required!")
	ENDIF()
ENDIF()

# Find bash itself ...
FIND_PROGRAM(BASH_TOOL_PATH bash
	HINTS ${GITHUB_BIN_DIR} DOC "The bash executable")

# Dumpbin is a windows dependency analaysis tool required for packaging
IF(WIN32 AND OGS_PACKAGING)
	FIND_PROGRAM(DUMPBIN_TOOL_PATH dumpbin DOC "Windows dependency analysis tool")
	IF(NOT DUMPBIN_TOOL_PATH)
		MESSAGE(FATAL_ERROR "Dumpbin was not found but is required for packaging!")
	ENDIF()
ENDIF()

########################
### Find other stuff ###
########################

# Check if on Jenkins
IF(NOT $ENV{JENKINS_URL} STREQUAL "")
	SET(JENKINS_URL $ENV{JENKINS_URL})
	SET(JENKINS_JOB_NAME $ENV{JOB_NAME})
ENDIF()

######################
### Find libraries ###
######################

# Do not search for libs if this option is set
IF(OGS_NO_EXTERNAL_LIBS)
	RETURN()
ENDIF() # OGS_NO_EXTERNAL_LIBS

# Clang does not have OpenMP support atm, see https://github.com/ufz/ogs/issues/8
IF(NOT COMPILER_IS_CLANG)
	FIND_PACKAGE(OpenMP)
ENDIF () # !clang
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

FIND_PACKAGE(Metis QUIET)

## Qt4 library ##
IF(NOT OGS_DONT_USE_QT)
	FIND_PACKAGE( Qt4 4.7)
ENDIF(NOT OGS_DONT_USE_QT)

IF ( QT4_FOUND )
	# Enable more modules
	SET(QT_USE_QTOPENGL TRUE)
	SET(QT_USE_QTSQL TRUE)
	SET(QT_USE_QTTEST TRUE)
	SET(QT_USE_QTXML TRUE)
	SET(QT_USE_QTXMLPATTERNS TRUE)
	INCLUDE( ${QT_USE_FILE} )
	ADD_DEFINITIONS(${QT_DEFINITIONS})
ENDIF (QT4_FOUND )

## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD ON )
FIND_PACKAGE ( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
	SET (HAVE_PTHREADS TRUE)
	ADD_DEFINITIONS(-DHAVE_PTHREADS)
ENDIF (CMAKE_USE_PTHREADS_INIT )

# blas
#FIND_PACKAGE ( BLAS QUIET )

# lapack
FIND_PACKAGE ( LAPACK QUIET )

## geotiff ##
IF(NOT MSVC)
	FIND_PACKAGE( LibTiff )
ENDIF() # NOT MSVC
FIND_PACKAGE( LibGeoTiff )
IF(libgeotiff_FOUND)
	ADD_DEFINITIONS(-Dlibgeotiff_FOUND)
ENDIF() # libgeotiff_FOUND

## lis ##
IF(OGS_USE_LIS)
	FIND_PACKAGE( LIS REQUIRED )
ENDIF()


IF(OGS_USE_PETSC)
    MESSAGE (STATUS  "Configuring for PETSc" )
   
    SET(OGS_USE_BOOSTMPI OFF)   
    SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake/findPETSC)
    FIND_PACKAGE(PETSc REQUIRED)
 
    include_directories(
              ${PETSC_INCLUDES} 
     )

    FIND_PACKAGE(MPI REQUIRED)

    ADD_DEFINITIONS(-DOGS_USE_PETSC)
ENDIF()
 