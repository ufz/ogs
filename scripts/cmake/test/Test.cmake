# Find tools and data
find_program(DIFF_TOOL_PATH diff)
find_program(TIME_TOOL_PATH time)
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
if(NOT GREP_TOOL_PATH)
    message(STATUS "grep-command is required for memcheck tester but was not found! All corresponding tests are disabled.")
endif()

enable_testing() # Enable CTest

# See http://www.vtk.org/Wiki/CMake/Testing_With_CTest for some customization options
set(CTEST_CUSTOM_TESTS_IGNORE test-harness) # ignore logog test
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/CTestCustom.cmake.in
    ${PROJECT_BINARY_DIR}/CTestCustom.cmake
)

include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/AddTest.cmake)

set(NUM_CTEST_PROCESSORS 3)

if(CMAKE_CONFIGURATION_TYPES)
    set(CONFIG_PARAMETER --build-config "$<CONFIGURATION>")
endif()
add_custom_target(ctest-cleanup ${CMAKE_COMMAND} -E remove Tests/ctest.log)
add_custom_target(
    ctest
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest.log
    --exclude-regex LARGE
    ${CONFIG_PARAMETER} --parallel ${NUM_CTEST_PROCESSORS} --test-action test
    DEPENDS ogs vtkdiff ctest-cleanup
)
add_custom_target(
    ctest-serial
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest.log
    --exclude-regex LARGE
    ${CONFIG_PARAMETER} --test-action test
    DEPENDS ogs vtkdiff ctest-cleanup
)
add_custom_target(ctest-large-cleanup ${CMAKE_COMMAND} -E remove Tests/ctest-large.log)
add_custom_target(
    ctest-large
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest-large.log
    ${CONFIG_PARAMETER} --parallel ${NUM_CTEST_PROCESSORS} --test-action test
    DEPENDS ogs vtkdiff ctest-large-cleanup
)
add_custom_target(
    ctest-large-serial
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest-large.log
    ${CONFIG_PARAMETER} --test-action test
    DEPENDS ogs vtkdiff ctest-large-cleanup
)
set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES ${PROJECT_BINARY_DIR}/Tests/lfs-data
)

set_target_properties(ctest PROPERTIES FOLDER Testing)
set_target_properties(ctest-large PROPERTIES FOLDER Testing)
set_target_properties(ctest-large-serial PROPERTIES FOLDER Testing)
set_target_properties(ctest-cleanup PROPERTIES FOLDER Testing)
set_target_properties(ctest-large-cleanup PROPERTIES FOLDER Testing)
