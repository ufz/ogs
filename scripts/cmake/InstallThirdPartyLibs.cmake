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
if(NOT INSTALL_PREFIX OR NOT GENERATOR)
	message(FATAL_ERROR "You need to specify an INSTALL_PREFIX and a GENERATOR: cmake -DINSTALL_PREFIX=C:/libs -DGENERATOR=\"Visual Studio 11 Win64\" -P ${CMAKE_CURRENT_LIST_DIR}/InstallThirdPartyLibs.cmake")
endif()
file(TO_CMAKE_PATH ${INSTALL_PREFIX} INSTALL_PREFIX)
file(MAKE_DIRECTORY ${INSTALL_PREFIX})
if(NOT IS_DIRECTORY ${INSTALL_PREFIX})
	message(FATAL_ERROR "Directory ${INSTALL_PREFIX} is not writable!")
endif()
if(NOT ${GENERATOR} STREQUAL "Visual Studio 11" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 11 Win64" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 12" AND
   NOT ${GENERATOR} STREQUAL "Visual Studio 12 Win64" AND
   NOT ${GENERATOR} STREQUAL "Unix Makefiles")
	message(FATAL_ERROR "Make sure to specify a supported GENERATOR!")
endif()

# CMake setup
set(CMAKE_PREFIX_PATH ${INSTALL_PREFIX})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
	"${CMAKE_CURRENT_LIST_DIR}/cmake"
	"${CMAKE_CURRENT_LIST_DIR}")
include(ThirdPartyLibVersions)

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(MAKE_PARALLEL_ARGS "")
if(NOT WIN32 AND NOT NUM_PROCESSORS EQUAL 0)
	set(MAKE_PARALLEL_ARGS "--" "-j${NUM_PROCESSORS}")
endif()

set(VISUAL_STUDIO_PARALLEL "")
if(WIN32)
	set(VISUAL_STUDIO_PARALLEL "-DCMAKE_CXX_FLAGS=\"/MP\"")
endif()

# Eigen
file(DOWNLOAD ${OGS_EIGEN_URL} tmp/eigen.tar.gz
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_EIGEN_MD5})
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf eigen.tar.gz
	WORKING_DIRECTORY tmp)
file(GLOB EIGEN_SOURCE_DIR tmp/eigen-*)
file(MAKE_DIRECTORY tmp/build-eigen)
execute_process(COMMAND ${CMAKE_COMMAND} ${EIGEN_SOURCE_DIR} -G "${GENERATOR}" -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
	WORKING_DIRECTORY tmp/build-eigen)
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install
	WORKING_DIRECTORY tmp/build-eigen)


# VTK
file(DOWNLOAD ${OGS_VTK_URL} tmp/vtk.tar.gz
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_VTK_MD5})
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf vtk.tar.gz
	WORKING_DIRECTORY tmp)
file(GLOB VTK_SOURCE_DIR tmp/VTK-*)
file(MAKE_DIRECTORY tmp/build-vtk)
execute_process(COMMAND ${CMAKE_COMMAND} ${VTK_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} -DBUILD_TESTING=OFF -DModule_vtkGUISupportQtOpenGL=ON
	-DBUILD_SHARED_LIBS=OFF -DCMAKE_DEBUG_POSTFIX=d ${VISUAL_STUDIO_PARALLEL}
	WORKING_DIRECTORY tmp/build-vtk)
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-vtk)
if(WIN32)
	execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Debug --target install
		WORKING_DIRECTORY tmp/build-vtk)
endif()

# GeoTiff
file(DOWNLOAD ${OGS_TIFF_URL} tmp/tiff.zip
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_TIFF_MD5})
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf tiff.zip
	WORKING_DIRECTORY tmp)
file(GLOB TIFF_SOURCE_DIR tmp/tiff-*)
file(MAKE_DIRECTORY tmp/build-tiff)
execute_process(COMMAND ${CMAKE_COMMAND} ${TIFF_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} ${VISUAL_STUDIO_PARALLEL}
	WORKING_DIRECTORY tmp/build-tiff)
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-tiff)

file(DOWNLOAD ${OGS_GEOTIFF_URL} tmp/geotiff.zip
	SHOW_PROGRESS EXPECTED_MD5 ${OGS_GEOTIFF_MD5})
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf geotiff.zip
	WORKING_DIRECTORY tmp)
file(GLOB GEOTIFF_SOURCE_DIR tmp/geotiff-*)
file(MAKE_DIRECTORY tmp/build-geotiff)
execute_process(COMMAND ${CMAKE_COMMAND} ${GEOTIFF_SOURCE_DIR} -G "${GENERATOR}"
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} -DWITH_PROJ4=OFF ${VISUAL_STUDIO_PARALLEL}
	WORKING_DIRECTORY tmp/build-geotiff)
execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release --target install ${MAKE_PARALLEL_ARGS}
	WORKING_DIRECTORY tmp/build-geotiff)

# Cleanup
file(REMOVE_RECURSE tmp)

message(STATUS "Finished!")
if(WIN32)
	file(TO_NATIVE_PATH ${INSTALL_PREFIX} INSTALL_PREFIX)
	message(STATUS "Now make sure to create an environment variable CMAKE_LIBRARY_SEARCH_PATH which points to ${INSTALL_PREFIX}!")
	message(STATUS "Make also sure to append %CMAKE_LIBRARY_SEARCH_PATH%\\bin to your PATH environment variable!")
endif()
