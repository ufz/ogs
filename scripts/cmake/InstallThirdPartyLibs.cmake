# - Install third-party libraries system-wide (intended for Windows)
#
# 2014-02-25, Lars Bilke
#
# Currently this will build:
#   Eigen, VTK, GeoTiff
#
# Usage:
#   Run this from wherever you want, it will create a tmp-directory, builds all
#   libraries in it and installs them at the given INSTALL_PREFIX location. The
#   tmp directory is deleted after finishing. Make sure to also pass your
#   GENERATOR.
#
#     cmake -DINSTALL_PREFIX=C:/libs -DGENERATOR="Visual Studio 12 Win64" \
#       -P path_to_ogs_source/scripts/cmake/InstallThirdPartyLibs.cmake
#
# Supported Generators:
#   - Visual Studio 11
#   - Visual Studio 11 Win64
#   - Visual Studio 12
#   - Visual Studio 12 Win64
#   - Unix Makefiles
#

# Argument checking
IF(NOT INSTALL_PREFIX OR NOT GENERATOR)
	MESSAGE(FATAL_ERROR "You need to specify an INSTALL_PREFIX and a GENERATOR: cmake -DINSTALL_PREFIX=C:/libs -DGENERATOR=\"Visual Studio 11 Win64\" -P ${CMAKE_CURRENT_LIST_DIR}/InstallThirdPartyLibs.cmake")
ENDIF()
FILE(TO_CMAKE_PATH ${INSTALL_PREFIX} INSTALL_PREFIX)
FILE(MAKE_DIRECTORY ${INSTALL_PREFIX})
IF(NOT IS_DIRECTORY ${INSTALL_PREFIX})
	MESSAGE(FATAL_ERROR "Directory ${INSTALL_PREFIX} is not writable!")
ENDIF()
IF(NOT ${GENERATOR} STREQUAL "Visual Studio 11" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 11 Win64" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 12" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 12 Win64" AND
   NOT ${GENERATOR} STREQUAL "Unix Makefiles")
	MESSAGE(FATAL_ERROR "Make sure to specify a supported GENERATOR!")
ENDIF()

# CMake setup
SET(CMAKE_PREFIX_PATH ${INSTALL_PREFIX})
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
	"${CMAKE_CURRENT_LIST_DIR}/cmake"
	"${CMAKE_CURRENT_LIST_DIR}")
INCLUDE(ThirdPartyLibVersions)

INCLUDE(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
SET(MAKE_PARALLEL_ARGS "")
IF(NOT WIN32 AND NOT NUM_PROCESSORS EQUAL 0)
	SET(MAKE_PARALLEL_ARGS "--" "-j${NUM_PROCESSORS}")
ENDIF()

# Eigen
FILE(DOWNLOAD ${OGS_EIGEN_URL} tmp/eigen.tar.gz
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_EIGEN_MD5})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xf eigen.tar.gz
	WORKING_DIRECTORY tmp)
FILE(GLOB EIGEN_SOURCE_DIR tmp/eigen-*)
FILE(MAKE_DIRECTORY tmp/build-eigen)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} ${EIGEN_SOURCE_DIR} -G "${GENERATOR}" -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
	WORKING_DIRECTORY tmp/build-eigen)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install
	WORKING_DIRECTORY tmp/build-eigen)


# VTK
FILE(DOWNLOAD ${OGS_VTK_URL} tmp/vtk.tar.gz
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_VTK_MD5})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xf vtk.tar.gz
	WORKING_DIRECTORY tmp)
FILE(GLOB VTK_SOURCE_DIR tmp/VTK-*)
FILE(MAKE_DIRECTORY tmp/build-vtk)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} ${VTK_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} -DBUILD_TESTING=OFF -DModule_vtkGUISupportQtOpenGL=ON
	-DBUILD_SHARED_LIBS=OFF -DCMAKE_DEBUG_POSTFIX=d
	WORKING_DIRECTORY tmp/build-vtk)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-vtk)
IF(WIN32)
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install
		WORKING_DIRECTORY tmp/build-vtk)
ENDIF()

# GeoTiff
FILE(DOWNLOAD ${OGS_TIFF_URL} tmp/tiff.zip
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_TIFF_MD5})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xf tiff.zip
	WORKING_DIRECTORY tmp)
FILE(GLOB TIFF_SOURCE_DIR tmp/tiff-*)
FILE(MAKE_DIRECTORY tmp/build-tiff)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} ${TIFF_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
	WORKING_DIRECTORY tmp/build-tiff)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-tiff)

FILE(DOWNLOAD ${OGS_GEOTIFF_URL} tmp/geotiff.zip
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_GEOTIFF_MD5})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xf geotiff.zip
	WORKING_DIRECTORY tmp)
FILE(GLOB GEOTIFF_SOURCE_DIR tmp/geotiff-*)
FILE(MAKE_DIRECTORY tmp/build-geotiff)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} ${GEOTIFF_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} -DWITH_PROJ4=OFF
	WORKING_DIRECTORY tmp/build-geotiff)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-geotiff)

# Cleanup
FILE(REMOVE_RECURSE tmp)

MESSAGE(STATUS "Finished!")
IF(WIN32)
	MESSAGE(STATUS "Now make sure to create an environment variable CMAKE_LIBRARY_SEARCH_PATH which points to ${INSTALL_PREFIX}!")
	MESSAGE(STATUS "Make also sure to append %CMAKE_LIBRARY_SEARCH_PATH%\\bin to your PATH environment variable!")
ENDIF()
