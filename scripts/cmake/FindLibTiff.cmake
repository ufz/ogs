# - Try to find libtiff
# Once done, this will define
#
#  libtiff_FOUND
#  libtiff_INCLUDE_DIRS
#  libtiff_LIBRARIES

if (NOT libtiff_FOUND)

	include(LibFindMacros)

	find_path( libtiff_INCLUDE_DIR
		NAMES tiff.h
		PATHS
			/usr/include
			${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
			$ENV{OGS_LIBS}/libtiff
			${OGS_LIBS_DIR}/libtiff/libtiff
		)

	if ( UNIX )
		find_library(libtiff_LIBRARIES
			NAMES tiff
			PATHS
				/usr/lib64
				/usr/lib
				${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
				${OGS_LIBS_DIR}/libtiff/libtiff
			)
	else ( UNIX )
		find_library(libtiff_LIBRARIES
			NAMES libtiff
			PATHS
				${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
				$ENV{OGS_LIBS}/libtiff
				${OGS_LIBS_DIR}/libtiff
			)
	endif ( UNIX )


	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT libtiff_LIBRARIES STREQUAL "libtiff_LIBRARIES-NOTFOUND" AND NOT libtiff_INCLUDE_DIR STREQUAL "libtiff_INCLUDE_DIR-NOTFOUND")
		set(libtiff_PROCESS_INCLUDES libtiff_INCLUDE_DIR)
		set(libtiff_PROCESS_LIBS libtiff_LIBRARIES)
		libfind_process(libtiff)
	else (NOT libtiff_LIBRARIES STREQUAL "libtiff_LIBRARIES-NOTFOUND" AND NOT libtiff_INCLUDE_DIR STREQUAL "libtiff_INCLUDE_DIR-NOTFOUND")
		message (STATUS "Could NOT find libtiff.")
	endif (NOT libtiff_LIBRARIES STREQUAL "libtiff_LIBRARIES-NOTFOUND" AND NOT libtiff_INCLUDE_DIR STREQUAL "libtiff_INCLUDE_DIR-NOTFOUND")

endif (NOT libtiff_FOUND)
