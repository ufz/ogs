# - Try to find libgeotiff
#
# Search directory
#
#  GEOTIFF_DIR
#
# Once done, this will define
#
#  GEOTIFF_FOUND
#  GEOTIFF_INCLUDE_DIRS
#  GEOTIFF_LIBRARIES

set(SEARCH_DIRS
	${GEOTIFF_DIR}
	$ENV{GEOTIFF_DIR}
	${CMAKE_SOURCE_DIR}/../Libs/libgeotiff
	$ENV{OGS_LIBS}/libgeotiff
	${OGS_LIBS_DIR}/libgeotiff
)

find_path( libgeotiff_INCLUDE_DIR
	NAMES geotiff.h
	PATHS
		/usr/include
		/usr/include/libgeotiff
		/usr/include/geotiff
		/usr/local/include
		/opt/boxen/homebrew/include
		${SEARCH_DIRS}
)

find_library(libgeotiff_LIBRARY
	NAMES geotiff
	PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib
		/opt/boxen/homebrew/lib
		${SEARCH_DIRS}
)

###
# Dependencies
###
set(_deps_libs)
set(_deps_includes)
set(_deps_check)

if(NOT MSVC)
	find_package(LibTiff)
endif()
list(APPEND _deps_libs ${TIFF_LIBRARIES})
list(APPEND _deps_includes ${TIFF_INCLUDE_DIRS})
list(APPEND _deps_check TIFF_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GEOTIFF
	REQUIRED_VARS
	libgeotiff_LIBRARY
	libgeotiff_INCLUDE_DIR
	${_deps_check}
)

if(GEOTIFF_FOUND)
	set(GEOTIFF_INCLUDE_DIRS ${libgeotiff_INCLUDE_DIR} ${_deps_includes})
	set(GEOTIFF_LIBRARIES ${libgeotiff_LIBRARY} ${_deps_libs})
endif()
