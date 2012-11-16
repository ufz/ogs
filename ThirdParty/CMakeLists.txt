# Add logog subdirectory and group its targets in a Visual Studio folder
ADD_SUBDIRECTORY(logog)
IF(BUILD_TESTING)
	SET_PROPERTY(TARGET Continuous PROPERTY FOLDER "logog")
ENDIF()
SET_PROPERTY(TARGET logog PROPERTY FOLDER "logog")
SET_PROPERTY(TARGET test-logog PROPERTY FOLDER "logog")

ADD_SUBDIRECTORY(zlib)

ADD_SUBDIRECTORY(gtest)

# This is a workaround for Travis builds.
OPTION(OGS_DONT_USE_OPTIONAL_SUBMODULES "Option for disabling optional submodules" OFF)

# VtkFbxConverter
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter/CMakeLists.txt
	AND NOT OGS_DONT_USE_OPTIONAL_SUBMODULES)
	ADD_SUBDIRECTORY(VtkFbxConverter)
ENDIF()

# VtkOsgConverter
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkOsgConverter/CMakeLists.txt
	AND NOT OGS_DONT_USE_OPTIONAL_SUBMODULES)
	ADD_SUBDIRECTORY(VtkOsgConverter)
ENDIF()