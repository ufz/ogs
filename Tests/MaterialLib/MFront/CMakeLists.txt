mfront_behaviours_check_library(
    OgsMFrontBehaviourForUnitTests
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
