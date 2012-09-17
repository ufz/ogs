# - Try to find Shapelib
# Once done, this will define
#
#  Shapelib_FOUND
#  Shapelib_INCLUDE_DIRS
#  Shapelib_LIBRARIES

if (NOT Shapelib_FOUND)

	include(LibFindMacros)

	set(SEARCH_DIRS
		${CMAKE_SOURCE_DIR}/../Libs/shapelib
		C:/OGS_Libs/shapelib
		$ENV{OGS_LIBS}/shapelib
		${OGS_LIBS_DIR_FOUND}/shapelib
	)

	find_path( Shapelib_INCLUDE_DIR
		NAMES shapefil.h
		PATHS
			/usr/include/libshp
			${SEARCH_DIRS}
		)

	find_library(Shapelib_LIBRARIES
		NAMES shp
		PATHS ${SEARCH_DIRS}
	)
	find_library(Shapelib_LIBRARIES
		NAMES shapelib
		PATHS ${SEARCH_DIRS}
	)

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT Shapelib_LIBRARIES STREQUAL "Shapelib_LIBRARIES-NOTFOUND" AND NOT Shapelib_INCLUDE_DIR STREQUAL "Shapelib_INCLUDE_DIR-NOTFOUND")
		set(Shapelib_PROCESS_INCLUDES Shapelib_INCLUDE_DIR)
		set(Shapelib_PROCESS_LIBS Shapelib_LIBRARIES)
		libfind_process(Shapelib)
	else (NOT Shapelib_LIBRARIES STREQUAL "Shapelib_LIBRARIES-NOTFOUND" AND NOT Shapelib_INCLUDE_DIR STREQUAL "Shapelib_INCLUDE_DIR-NOTFOUND")
		#message (STATUS "Warning: shapelib not found!")
	endif (NOT Shapelib_LIBRARIES STREQUAL "Shapelib_LIBRARIES-NOTFOUND" AND NOT Shapelib_INCLUDE_DIR STREQUAL "Shapelib_INCLUDE_DIR-NOTFOUND")

endif (NOT Shapelib_FOUND)
