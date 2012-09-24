IF(OGS_COVERAGE)
	SETUP_TARGET_FOR_COVERAGE_COBERTURA(ctest_coverage ctest "ctest_coverage_results" "-j;${PROCESSOR_COUNT}")
ENDIF()