# - Try to find libgeotiff
#
# Search directory
#
#  libgeotiff_DIR
#
# Once done, this will define
#
#  libgeotiff_FOUND
#  libgeotiff_INCLUDE_DIRS
#  libgeotiff_LIBRARIES

if (NOT libgeotiff_FOUND)

	include(LibFindMacros)

	set(SEARCH_DIRS
		${libgeotiff_DIR}
		$ENV{libgeotiff_DIR}
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
				${SEARCH_DIRS}
			)

	find_library(libgeotiff_LIBRARIES
		NAMES geotiff
		PATHS
			/usr/lib64
			/usr/lib
			${SEARCH_DIRS}
		)


	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT libgeotiff_LIBRARIES STREQUAL "libgeotiff_LIBRARIES-NOTFOUND" AND NOT libgeotiff_INCLUDE_DIR STREQUAL "libgeotiff_INCLUDE_DIR-NOTFOUND")
		set(libgeotiff_PROCESS_INCLUDES libgeotiff_INCLUDE_DIR)
		set(libgeotiff_PROCESS_LIBS libgeotiff_LIBRARIES)
		libfind_process(libgeotiff)
	else (NOT libgeotiff_LIBRARIES STREQUAL "libgeotiff_LIBRARIES-NOTFOUND" AND NOT libgeotiff_INCLUDE_DIR STREQUAL "libgeotiff_INCLUDE_DIR-NOTFOUND")
		message (STATUS "Could NOT find libgeotiff.")
	endif (NOT libgeotiff_LIBRARIES STREQUAL "libgeotiff_LIBRARIES-NOTFOUND" AND NOT libgeotiff_INCLUDE_DIR STREQUAL "libgeotiff_INCLUDE_DIR-NOTFOUND")

	SET ( libgeotiff_INCLUDE_DIR ${libgeotiff_INCLUDE_DIR} ${libtiff_INCLUDE_DIR} CACHE STRING "libgeotiff include directories." FORCE )
	SET ( libgeotiff_LIBRARIES ${libgeotiff_LIBRARIES} ${libtiff_LIBRARIES} CACHE STRING "libgeotiff link libraries." FORCE )

endif (NOT libgeotiff_FOUND)
