# Add logog subdirectory and group its targets in a Visual Studio folder
ADD_SUBDIRECTORY(logog)
IF(BUILD_TESTING)
	SET_PROPERTY(TARGET Continuous PROPERTY FOLDER "logog")
ENDIF()
SET_PROPERTY(TARGET logog PROPERTY FOLDER "logog")
SET_PROPERTY(TARGET test-logog PROPERTY FOLDER "logog")

ADD_SUBDIRECTORY(zlib)

ADD_SUBDIRECTORY(gtest)

# VtkFbxConverter
#SET(VtkFbxConverter_DIR ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter)
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter/CMakeLists.txt)
	ADD_SUBDIRECTORY(VtkFbxConverter)
	SET(VTKFBXCONVERTER_FOUND ON CACHE BOOL "Is the VtkFbxConverter submodule checked out?")
ENDIF()