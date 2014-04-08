# Find tools and data
FIND_PROGRAM(DIFF_TOOL_PATH diff)
FIND_PROGRAM(GREP_TOOL_PATH grep)
FIND_PROGRAM(BASH_TOOL_PATH bash)
FIND_PROGRAM(VALGRIND_TOOL_PATH valgrind)
FIND_FILE(OGS-DATA_PATH ".ogs6-data.dummy"
	HINTS ${CMAKE_SOURCE_DIR}/../ogs6-data)
GET_FILENAME_COMPONENT(OGS-DATA_PATH ${OGS-DATA_PATH} PATH)

ENABLE_TESTING() # Enable CTest

# See http://www.vtk.org/Wiki/CMake/Testing_With_CTest for some customization options
SET(CTEST_CUSTOM_TESTS_IGNORE test-harness) # ignore logog test
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/CTestCustom.cmake.in
	${CMAKE_BINARY_DIR}/CTestCustom.cmake
)
