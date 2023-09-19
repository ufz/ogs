set(CTEST_PROJECT_NAME "ogs")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
set(CTEST_SUBMIT_URL
    "https://cdash.opengeosys.org/submit.php?project=${CTEST_PROJECT_NAME}"
)
set(CTEST_SUBMIT_INACTIVITY_TIMEOUT 30)
if(DEFINED ENV{CI_JOB_NAME})
    # Bug in CDash: Remove ":", see
    # https://github.com/Kitware/CDash/issues/1292
    string(REPLACE ":" "" _build_name "$ENV{CI_JOB_NAME}")
    set(BUILDNAME ${_build_name})
endif()
