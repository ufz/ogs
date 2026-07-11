mfront_behaviours_check_library(
    OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly
    BentoniteBehaviour
)

# Disable warnings for generated OgsMFrontBehaviourForUnitTests
target_compile_options(
    OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly
    PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
            $<$<CXX_COMPILER_ID:MSVC>:/W0>
)

# generalmod is imported third-party code (TRIAX, D. Masin); do not lint it.
set_target_properties(
    OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly PROPERTIES CXX_CLANG_TIDY ""
)

target_include_directories(
    OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

set_target_properties(OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly PROPERTIES EXCLUDE_FROM_ALL FALSE)
install(TARGETS OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly)
