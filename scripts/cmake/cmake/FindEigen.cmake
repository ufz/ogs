FIND_PATH(EIGEN_INCLUDE_DIR NAMES Eigen/Core
	HINTS ${PC_EIGEN_INCLUDEDIR} ${PC_EIGEN_INCLUDE_DIRS} "$ENV{PROGRAMFILES}/Eigen 3.0.0/include/eigen3"
	PATH_PREFIX eigen3
)

set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR)


mark_as_advanced(EIGEN_INCLUDE_DIR)

if(EIGEN_FOUND)
	INCLUDE_DIRECTORIES( ${EIGEN_INCLUDE_DIRS} )
	message(STATUS "Eigen found (include: ${EIGEN_INCLUDE_DIRS})")
endif(EIGEN_FOUND)
