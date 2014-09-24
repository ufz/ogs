INCLUDE(CodeCoverage)

SET(COVERAGE_EXCLUDES
	'/gpfs0/*'
	'/usr/*'
	'${CMAKE_BINARY_DIR}/*'
	'${CMAKE_SOURCE_DIR}/Tests/*'
	'${CMAKE_SOURCE_DIR}/ThirdParty/*'
)

SETUP_TARGET_FOR_COVERAGE_COBERTURA(testrunner_coverage_cobertura testrunner "testrunner_coverage_cobertura_results" "-j;${PROCESSOR_COUNT}")
SETUP_TARGET_FOR_COVERAGE(testrunner_coverage testrunner "testrunner_coverage_results" "-j;${PROCESSOR_COUNT}")
IF(OGS_BUILD_GUI)
	SETUP_TARGET_FOR_COVERAGE(DataExplorer_coverage DataExplorer "DataExplorer_coverage_results")
ENDIF()
