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

find_package(Git REQUIRED)
set(GIT_TOOL_PATH ${GIT_EXECUTABLE} CACHE FILEPATH "The git command line interface" FORCE)

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

## CCache
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
	set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
	if(COMPILER_IS_CLANG)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
	endif()
endif(CCACHE_FOUND)

######################
### Find libraries ###
######################

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
	find_package( Qt5 5.5 REQUIRED Gui Widgets Xml XmlPatterns)

	cmake_policy(SET CMP0020 NEW)
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

if(OGS_USE_PETSC)
	message(STATUS "Configuring for PETSc")

	option(FORCE_PETSC_EXECUTABLE_RUNS "Force CMake to accept a given PETSc configuration" ON)

	##Force CMake to accept a given PETSc configuration in case the failure of MPI tests
	##This may cause the compilation broken.
	if(FORCE_PETSC_EXECUTABLE_RUNS)
		set(PETSC_EXECUTABLE_RUNS YES)
	endif()

	find_package(PETSc REQUIRED)

	include_directories(SYSTEM ${PETSC_INCLUDES})

	add_definitions(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)

endif()

## Check MPI package
if(OGS_USE_MPI)
	find_package(MPI REQUIRED)
	include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
endif()
