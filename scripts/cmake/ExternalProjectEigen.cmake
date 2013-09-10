MESSAGE(STATUS "Eigen inc: ${EIGEN3_INCLUDE_DIR}")

IF(EIGEN3_FOUND)
	RETURN()
ENDIF()

# First check for system Eigen
IF(NOT EIGEN3_INCLUDE_DIR)
	FIND_PACKAGE(Eigen3)
	IF(EIGEN3_FOUND)
		SET(EIGEN3_FOUND TRUE CACHE BOOL "Was Eigen found?" FORCE)
		SET(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR}" CACHE STRING "Eigen include dir" FORCE)
		RETURN()
	ELSE()
		SET(EIGEN3_INCLUDE_DIR "")
	ENDIF()
ENDIF()

ExternalProject_Add(Eigen
	PREFIX ${CMAKE_BINARY_DIR}/External/eigen
	URL http://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz
	URL_MD5 9559c34af203dde5f3f1d976d859c5b3
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)
ExternalProject_Get_Property( Eigen source_dir )

IF(NOT EIGEN3_INCLUDE_DIR)
	SET(EIGEN3_FOUND TRUE CACHE BOOL "Was Eigen found?" FORCE)
	SET( EIGEN3_INCLUDE_DIR ${source_dir} CACHE INTERNAL "Eigen include dir")
	MESSAGE(STATUS "Downloading Eigen automatically.")
ENDIF()
