# Disable all warnings
if(COMPILER_IS_CLANG OR COMPILER_IS_GCC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
endif()

# Add logog subdirectory and group its targets in a Visual Studio folder
add_subdirectory(logog)
if(BUILD_TESTING)
	set_property(TARGET Continuous PROPERTY FOLDER "logog")
endif()
set_property(TARGET logog PROPERTY FOLDER "logog")
set_property(TARGET test-logog PROPERTY FOLDER "logog")
set_target_properties(test-logog PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_subdirectory(zlib)

add_subdirectory(gtest)

add_subdirectory(shape)

# This is a workaround for Travis builds.
option(OGS_USE_OPTIONAL_SUBMODULES "Option for enabling optional submodules" OFF)

# VtkFbxConverter
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter/CMakeLists.txt
	AND OGS_USE_OPTIONAL_SUBMODULES)
	add_subdirectory(VtkFbxConverter)
endif()

# VtkOsgConverter
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkOsgConverter/CMakeLists.txt
	AND OGS_USE_OPTIONAL_SUBMODULES)
	add_subdirectory(VtkOsgConverter)
endif()

# tetgen
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tetgen/CMakeLists.txt)
	add_subdirectory(tetgen)
	if(APPLE AND OGS_BUILD_GUI)
		install (PROGRAMS $<TARGET_FILE:tetgen> DESTINATION DataExplorer.app/Contents/MacOS COMPONENT ogs_extras)
	else()
		install (PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin COMPONENT ogs_extras)
	endif()
endif()

# vtkdiff
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/vtkdiff/CMakeLists.txt)
	add_subdirectory(vtkdiff)
	install (PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin COMPONENT ogs_extras)
endif()
