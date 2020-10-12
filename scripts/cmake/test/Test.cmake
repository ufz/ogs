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

include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/AddTest.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/MeshTest.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/test/OgsTest.cmake)

if(CMAKE_CONFIGURATION_TYPES)
    set(CONFIG_PARAMETER --build-config "$<CONFIGURATION>")
endif()
add_custom_target(ctest-cleanup ${CMAKE_COMMAND} -E remove -f Tests/ctest.log)

set(test_dependencies ogs vtkdiff)
if(OGS_BUILD_UTILS)
    list(APPEND test_dependencies
        partmesh
        MapGeometryToMeshSurface
        generateStructuredMesh
    )
endif()
if(OGS_USE_MFRONT)
    list(APPEND test_dependencies
        MFrontGenericBehaviourInterfaceTest
        MFrontGenericBehaviourInterfaceTest2
        MFrontGenericBehaviourInterfaceTest3
        BoundsCheckTest
        ParameterTest
        IntegrateTest
        IntegrateTest2
        IntegrateTest2b
        IntegrateTest3
        IntegrateTest3b
        BehaviourTest)
endif()

add_custom_target(
    ctest
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest.log
    --exclude-regex LARGE
    ${CONFIG_PARAMETER}
    --timeout 900 # 15 minutes
    DEPENDS ${test_dependencies} ctest-cleanup
    USES_TERMINAL
)

add_custom_target(ctest-large-cleanup ${CMAKE_COMMAND} -E remove -f Tests/ctest-large.log)

add_custom_target(
    ctest-large
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
    --force-new-ctest-process
    --output-on-failure --output-log Tests/ctest-large.log
    --tests-regex LARGE
    ${CONFIG_PARAMETER}
    --timeout 5400 # 90 minutes
    DEPENDS ${test_dependencies} ctest-large-cleanup
    USES_TERMINAL
)

set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES ${PROJECT_BINARY_DIR}/Tests/Data
)

set_target_properties(ctest ctest-large ctest-cleanup ctest-large-cleanup
    PROPERTIES FOLDER Testing)

add_dependencies(ctest ogs)
add_dependencies(ctest-large ogs)

configure_file(${PROJECT_SOURCE_DIR}/scripts/test/buildinfo.in.yaml ${PROJECT_BINARY_DIR}/buildinfo.yaml)
