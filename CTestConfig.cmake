set(CTEST_PROJECT_NAME "ogs")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
set(CTEST_SUBMIT_URL
    "https://cdash.opengeosys.org/submit.php?project=${CTEST_PROJECT_NAME}"
)
if(DEFINED ENV{CI_JOB_NAME})
    set(BUILDNAME $ENV{CI_JOB_NAME})
endif()