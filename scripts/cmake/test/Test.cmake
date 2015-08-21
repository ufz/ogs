# Find tools and data
find_program(DIFF_TOOL_PATH diff)
find_program(NUMDIFF_TOOL_PATH numdiff)
# find_program(TIME_TOOL_PATH time) # TODO: does not work Travis
set(TIME_TOOL_PATH time)
find_program(GREP_TOOL_PATH grep)
find_program(BASH_TOOL_PATH bash)
find_program(VALGRIND_TOOL_PATH valgrind)
find_program(MPIRUN_TOOL_PATH mpirun)

if(NOT TIME_TOOL_PATH)
	message(STATUS "time-command is required for time wrapper but was not found! All corresponding tests are disabled.")
endif()
if(NOT VALGRIND_TOOL_PATH)
	message(STATUS "Valgrind is required for memcheck wrapper but was not found! All corresponding tests are disabled.")
endif()
if(NOT VALGRIND_TOOL_PATH)
	message(STATUS "Valgrind is required for callgrind wrapper but was not found! All corresponding tests are disabled.")
endif()
if(NOT MPIRUN_TOOL_PATH)
	message(STATUS "mpirun is required for mpirun wrapper but was not found! All corresponding tests are disabled.")
endif()
if(NOT DIFF_TOOL_PATH)
	message(STATUS "diff-command is required for diff tester but was not found! All corresponding tests are disabled.")
endif()
if(NOT NUMDIFF_TOOL_PATH)
	message(STATUS "numdiff-command is required for numdiff tester but was not found! All corresponding tests are disabled.")
endif()
if(NOT GREP_TOOL_PATH)
	message(STATUS "grep-command is required for memcheck tester but was not found! All corresponding tests are disabled.")
endif()

enable_testing() # Enable CTest

# See http://www.vtk.org/Wiki/CMake/Testing_With_CTest for some customization options
set(CTEST_CUSTOM_TESTS_IGNORE test-harness) # ignore logog test
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/CTestCustom.cmake.in
	${CMAKE_BINARY_DIR}/CTestCustom.cmake
)

include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/AddTest.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/Data.cmake)

if(CMAKE_CONFIGURATION_TYPES)
	set(CONFIG_PARAMETER --build-config "$<CONFIGURATION>")
endif()
add_custom_target(
	ctest -T Test
	COMMAND ${CMAKE_CTEST_COMMAND}
	--force-new-ctest-process --output-on-failure --exclude-regex LARGE
	${CONFIG_PARAMETER} --parallel ${NUM_PROCESSORS} --test-action test
	DEPENDS data ogs vtkdiff
)
add_custom_target(
	ctest-large -T Test
	COMMAND ${CMAKE_CTEST_COMMAND}
	--force-new-ctest-process --output-on-failure --tests-regex LARGE
	${CONFIG_PARAMETER} --parallel ${NUM_PROCESSORS} --test-action test
	DEPENDS data ogs vtkdiff
)
