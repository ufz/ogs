# - Try to find libgeotiff
#
# Once done, this will define
#
#  GEOTIFF_FOUND
#  GEOTIFF_INCLUDE_DIRS
#  GEOTIFF_LIBRARIES

find_path( libgeotiff_INCLUDE_DIR geotiff.h)

find_library(libgeotiff_LIBRARY geotiff)

###
# Dependencies
###
set(_deps_libs)
set(_deps_includes)
set(_deps_check)

if(NOT MSVC)
	find_package(TIFF)
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
