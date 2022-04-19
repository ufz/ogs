include(FeatureSummary)
add_feature_info(OGS OGS_BUILD_CLI "The OGS simulator (OGS_BUILD_CLI)")
add_feature_info(
    Utilities OGS_BUILD_UTILS "Command line tools (OGS_BUILD_UTILS)"
)
add_feature_info(
    DataExplorer OGS_BUILD_GUI
    "Graphical data exploration and processing (OGS_BUILD_GUI)"
)
add_feature_info(
    MFront OGS_USE_MFRONT "MFront material models (OGS_USE_MFRONT)"
)
add_feature_info(
    Python-interface OGS_USE_PYTHON
    "Python boundary conditions and source terms (OGS_USE_PYTHON)"
)
add_feature_info(
    Python-environment
    OGS_USE_PIP
    "Automatic Python dependency handling with a virtual environment (OGS_USE_PIP)"
)
add_feature_info(PETSc OGS_USE_PETSC "Parallel processing (OGS_USE_PETSC)")
add_feature_info(
    Tests OGS_BUILD_TESTING "Unit and benchmarks tests (OGS_BUILD_TESTING)"
)
add_feature_info(
    build-shared BUILD_SHARED_LIBS "Shared libraries (BUILD_SHARED_LIBS)"
)
add_feature_info(
    build-unity OGS_USE_UNITY_BUILDS "Unity build (OGS_USE_UNITY_BUILDS)"
)

feature_summary(WHAT PACKAGES_FOUND ENABLED_FEATURES)
