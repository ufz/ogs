# - Try to find libtiff
# Once done, this will define
#
#  libtiff_FOUND
#  libtiff_INCLUDE_DIRS
#  libtiff_LIBRARIES

find_path( libtiff_INCLUDE_DIR
	NAMES tiff.h
	PATHS
		/usr/include
		/usr/local/include
		/opt/boxen/homebrew/include
		${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
		$ENV{OGS_LIBS}/libtiff
		${OGS_LIBS_DIR}/libtiff/libtiff
	)

if ( UNIX )
	find_library(libtiff_LIBRARY
		NAMES tiff
		PATHS
			/usr/lib64
			/usr/lib
			/usr/local/lib
			/opt/boxen/homebrew/lib
			${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
			${OGS_LIBS_DIR}/libtiff/libtiff
	)
else ( UNIX )
	find_library(libtiff_LIBRARY
		NAMES libtiff
		PATHS
			${CMAKE_SOURCE_DIR}/../Libs/libtiff/libtiff
			$ENV{OGS_LIBS}/libtiff
			${OGS_LIBS_DIR}/libtiff
		)
endif ( UNIX )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TIFF
	REQUIRED_VARS
	libtiff_LIBRARY
	libtiff_INCLUDE_DIR
	${_deps_check}
)

if(TIFF_FOUND)
	set(TIFF_INCLUDE_DIRS ${libtiff_INCLUDE_DIR})
	set(TIFF_LIBRARIES ${libtiff_LIBRARY})
endif()
