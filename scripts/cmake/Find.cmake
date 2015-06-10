# Add custom library install prefixes
list(APPEND CMAKE_PREFIX_PATH
	$ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
	$ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
	${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

######################
### Find tools     ###
######################

# Find doxygen
if(WIN32)
	find_program(DOXYGEN_DOT_EXECUTABLE NAMES dot PATHS "$ENV{ProgramFiles}/Graphviz*/bin")
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

find_package(GitHub)

find_program(GIT_TOOL_PATH git HINTS ${GITHUB_BIN_DIR} DOC "The git command line interface")
if(NOT GIT_TOOL_PATH)
	if(WIN32)
		message(FATAL_ERROR "Git not found! Please install GitHub for Windows or Git!")
	else()
		message(FATAL_ERROR "Git not found but is required!")
	endif()
else()
	set(GIT_EXECUTABLE ${GIT_TOOL_PATH} CACHE FILE_PATH "" FORCE)
	set(GIT_FOUND TRUE CACHE BOOL "" FORCE)
endif()

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

######################
### Find libraries ###
######################

# Do not search for libs if this option is set
if(OGS_NO_EXTERNAL_LIBS)
	return()
endif() # OGS_NO_EXTERNAL_LIBS

find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Metis QUIET)

## Qt4 library ##
if(NOT OGS_DONT_USE_QT)
	find_package( Qt4 4.7)
endif()

if(QT4_FOUND)
	# Enable more modules
	set(QT_USE_QTOPENGL TRUE)
	set(QT_USE_QTSQL TRUE)
	set(QT_USE_QTTEST TRUE)
	set(QT_USE_QTXML TRUE)
	set(QT_USE_QTXMLPATTERNS TRUE)
	include( ${QT_USE_FILE} )
	add_definitions(${QT_DEFINITIONS} -DQT4_FOUND)
endif()

## pthread ##
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
	set(HAVE_PTHREADS TRUE)
	add_definitions(-DHAVE_PTHREADS)
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

if(OGS_USE_PETSC)
	message(STATUS "Configuring for PETSc")

	option(FORCE_PETSC_EXECUTABLE_RUNS "Force CMake to accept a given PETSc configuration" ON)

	##Force CMake to accept a given PETSc configuration in case the failure of MPI tests
	##This may cause the compilation broken.
	if(FORCE_PETSC_EXECUTABLE_RUNS)
		set(PETSC_EXECUTABLE_RUNS YES)
	endif()

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/scripts/cmake/findPETSC")
	find_package(PETSc REQUIRED)

	include_directories(SYSTEM ${PETSC_INCLUDES})

	add_definitions(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)

endif()

## Check MPI package
if(OGS_USE_MPI)
	find_package(MPI REQUIRED)
	include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
endif()

