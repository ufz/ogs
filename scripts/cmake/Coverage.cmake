find_program(FASTCOV_PATH NAMES fastcov fastcov.py)
if(NOT FASTCOV_PATH AND NOT OGS_USE_PIP)
    message(
        FATAL_ERROR "Code coverage requires either fastcov or OGS_USE_PIP=ON."
    )
endif()

# https://github.com/linux-test-project/lcov/pull/125
if(APPLE)
    file(
        DOWNLOAD
        https://raw.githubusercontent.com/linux-test-project/lcov/41d8655951d6898511f98be2a2dbcfbe662f0b17/bin/genhtml
        ${PROJECT_BINARY_DIR}/bin/genhtml
    )
    set(GENHTML_PATH ${PROJECT_BINARY_DIR}/bin/genhtml)
endif()

include(CodeCoverage)
append_coverage_compiler_flags()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og")

# cmake-lint: disable=C0103
if(NOT FASTCOV_PATH)
    list(APPEND OGS_PYTHON_PACKAGES "fastcov==1.14")
    set(FASTCOV_PATH ${LOCAL_VIRTUALENV_BIN_DIR}/fastcov CACHE INTERNAL "")
endif()

if(DEFINED ENV{CI})
    set(COVERAGE_ADDITIONAL_ARGS SKIP_HTML)
endif()

# ~~~
# TODO: segfault in MeshLibMappedPropertyVector.Double|Int
# TODO: segfault in TestVtkMeshConverter.Conversion
# ~~~
setup_target_for_coverage_fastcov(
    NAME
    testrunner_coverage
    BASE_DIRECTORY
    ${PROJECT_BINARY_DIR}
    EXECUTABLE
    $<TARGET_FILE:testrunner>
    -l
    warn
    --gtest_filter=-MeshLibMappedPropertyVector.*:GeoLib.SearchNearestPointsInDenseGrid:TestVtkMeshConverter.Conversion
    DEPENDENCIES
    testrunner
    FASTCOV_ARGS
    --branch-coverage
    --include
    ${PROJECT_SOURCE_DIR}
    ${COVERAGE_ADDITIONAL_ARGS}
    EXCLUDE
    Applications/CLI/
    ProcessLib/
    Tests/
)

# TODO: segfault in Vtu2Grid
setup_target_for_coverage_fastcov(
    NAME
    ctest_coverage
    BASE_DIRECTORY
    ${PROJECT_BINARY_DIR}
    EXECUTABLE
    ctest
    -E
    "Vtu2Grid"
    DEPENDENCIES
    all
    FASTCOV_ARGS
    --branch-coverage
    --include
    ${PROJECT_SOURCE_DIR}
    ${COVERAGE_ADDITIONAL_ARGS}
    EXCLUDE
    Applications/CLI/
    Tests/
    POST_CMD
    perl
    -i
    -pe
    s!${PROJECT_SOURCE_DIR}/!!g
    ctest_coverage.json
    NO_DEMANGLE
)

if(UNIX)
    add_custom_target(clean_coverage find . -name '*.gcda' -delete)
endif()

configure_file(
    ${PROJECT_SOURCE_DIR}/scripts/test/generate_coverage_vis_data.in.py
    ${PROJECT_BINARY_DIR}/generate_coverage_vis_data.py @ONLY
)
