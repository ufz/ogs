# Find tools and data
FIND_PROGRAM(DIFF_TOOL_PATH diff)
# FIND_PROGRAM(TIME_TOOL_PATH time) # TODO: does not work Travis
SET(TIME_TOOL_PATH time)
FIND_PROGRAM(GREP_TOOL_PATH grep)
FIND_PROGRAM(BASH_TOOL_PATH bash)
FIND_PROGRAM(VALGRIND_TOOL_PATH valgrind)

ENABLE_TESTING() # Enable CTest

# See http://www.vtk.org/Wiki/CMake/Testing_With_CTest for some customization options
SET(CTEST_CUSTOM_TESTS_IGNORE test-harness) # ignore logog test
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/CTestCustom.cmake.in
	${CMAKE_BINARY_DIR}/CTestCustom.cmake
)

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/AddTest.cmake)
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/Data.cmake)

IF(CMAKE_CONFIGURATION_TYPES)
	ADD_CUSTOM_TARGET(
		ctest
		COMMAND ${CMAKE_CTEST_COMMAND}
		--force-new-ctest-process --output-on-failure
		--build-config "$<CONFIGURATION>"
		DEPENDS data
	)
ELSE()
	ADD_CUSTOM_TARGET(
		ctest
		COMMAND ${CMAKE_CTEST_COMMAND}
		--force-new-ctest-process --output-on-failure
		DEPENDS data
	)
ENDIF()
