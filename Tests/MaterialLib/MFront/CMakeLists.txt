mfront_behaviours_check_library(
    OgsMFrontBehaviourForUnitTests
    CheckIntegrationFailure
    CheckParamPassing1
    CheckStiffnessMatrixAllBlocks
    CheckStiffnessMatrixAllBlocksReordered
    CheckStiffnessMatrixExcessBlocks
    CheckStiffnessMatrixMissingBlocks
)

# Disable warnings for generated OgsMFrontBehaviourForUnitTests
target_compile_options(
    OgsMFrontBehaviourForUnitTests
    PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
            $<$<CXX_COMPILER_ID:MSVC>:/W0>
)

# Python generates the cases; CMake only registers the individual MTests.
set(_lode_generator
    "${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombLodeSmoothing/generate_mtests.py"
)
set_property(
    DIRECTORY
    APPEND
    PROPERTY
        CMAKE_CONFIGURE_DEPENDS
        "${_lode_generator}"
        "${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombLodeSmoothing/MohrCoulombAbboSloan_lode.mtest.in"
)
execute_process(
    COMMAND
        "${Python_EXECUTABLE}" "${_lode_generator}" --library
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libOgsMFrontBehaviour${CMAKE_SHARED_LIBRARY_SUFFIX}"
        --output-directory "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_VARIABLE _lode_tests
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY
)
foreach(lode_test IN LISTS _lode_tests)
    # TFEL may be built later as an external project; resolve mtest at test
    # time.
    add_test(
        NAME ${lode_test}
        COMMAND ${CMAKE_COMMAND} -E env "PATH=${TFELHOME}/bin:$ENV{PATH}" mtest
                --verbose=quiet --result-file-output=false "${lode_test}.mtest"
    )
    set_tests_properties(${lode_test} PROPERTIES LABELS "default;MFront;small")
endforeach()
