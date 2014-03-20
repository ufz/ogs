# Add custom library install prefixes
LIST(APPEND CMAKE_PREFIX_PATH
	$ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
	$ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
	${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

######################
### Find tools     ###
######################

# Find doxygen
IF(WIN32)
	FIND_PROGRAM(DOXYGEN_DOT_EXECUTABLE NAMES dot PATHS "$ENV{ProgramFiles}/Graphviz*/bin")
	FIND_PACKAGE(Doxygen QUIET)
	IF(DOXYGEN_DOT_PATH)
		FILE(TO_NATIVE_PATH ${DOXYGEN_DOT_PATH} DOXYGEN_DOT_PATH)
		SET(DOXYGEN_DOT_PATH "\"${DOXYGEN_DOT_PATH}\"")
	ENDIF()
ELSE()
	FIND_PACKAGE(Doxygen QUIET)
ENDIF()

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

# Dumpbin is a windows dependency analaysis tool required for packaging.
# Variable has to be named gp_cmd to override the outdated find routines
# of the GetPrerequisites CMake-module.
IF(WIN32)
	INCLUDE(MSVCPaths)
	FIND_PROGRAM(gp_cmd dumpbin DOC "Windows dependency analysis tool"
		PATHS ${MSVC_INSTALL_PATHS} PATH_SUFFIXES VC/bin)
	IF(gp_cmd)
		GET_FILENAME_COMPONENT(dir ${gp_cmd} PATH)
		SET(ENV{PATH} "${dir}/../../../Common7/IDE;$ENV{PATH}")
	ENDIF()
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
FIND_PACKAGE( LibGeoTiff )
IF(GEOTIFF_FOUND)
	ADD_DEFINITIONS(-DGEOTIFF_FOUND)
ENDIF() # GEOTIFF_FOUND

## lis ##
IF(OGS_USE_LIS)
	FIND_PACKAGE( LIS REQUIRED )
ENDIF()


IF(OGS_USE_PETSC)
    MESSAGE (STATUS  "Configuring for PETSc" )
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/scripts/cmake/findPETSC")
    FIND_PACKAGE(PETSc REQUIRED)

    include_directories(
              ${PETSC_INCLUDES}
     )

    FIND_PACKAGE(MPI REQUIRED)
ENDIF()
 
