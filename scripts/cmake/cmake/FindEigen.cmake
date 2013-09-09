# - Find Eigen3 includes
# Once done, this will define
#
#  EIGEN_INCLUDE_DIR - Eigen include directory

find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core
	PATH_PREFIX eigen3
	PATHS
		/usr/local/include
		/usr/include
	)

set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR)

# Mark the named cached variables as advanced
mark_as_advanced(EIGEN_INCLUDE_DIR)

if(EIGEN_FOUND)
	INCLUDE_DIRECTORIES( ${EIGEN_INCLUDE_DIRS} )
	message(STATUS "Eigen found (include: ${EIGEN_INCLUDE_DIRS})")
endif(EIGEN_FOUND)
