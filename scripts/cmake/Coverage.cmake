find_program(FASTCOV_PATH NAMES fastcov fastcov.py)
if(NOT FASTCOV_PATH AND NOT POETRY)
    message(FATAL_ERROR "Code coverage requires either fastcov or poetry.")
endif()

include(CodeCoverage)
append_coverage_compiler_flags()
# TODO: segfault in MeshLibMappedPropertyVector.Double?
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og")

if(NOT FASTCOV_PATH)
    execute_process(COMMAND ${CMD_COMMAND} poetry add fastcov==1.10
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    find_program(FASTCOV_PATH NAMES fastcov
        HINTS ${LOCAL_VIRTUALENV_BIN_DIRS}
        REQUIRED NO_DEFAULT_PATH
    )
endif()
setup_target_for_coverage_fastcov(
    NAME testrunner_coverage
    EXECUTABLE $<TARGET_FILE:testrunner> -l warn --gtest_filter=-MeshLibMappedPropertyVector.*:GeoLib.SearchNearestPointsInDenseGrid
    DEPENDENCIES testrunner
    FASTCOV_ARGS --include ${PROJECT_SOURCE_DIR}
    EXCLUDE
        Applications/CLI/
        ProcessLib/
        Tests/
)
