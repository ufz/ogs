INCLUDE(CodeCoverage)

SET(COVERAGE_EXCLUDES
	'/usr/*'
	'${CMAKE_BINARY_DIR}/*'
	'${CMAKE_SOURCE_DIR}/tests/*'
	'${CMAKE_SOURCE_DIR}/ThirdParty/*'
)

IF(JENKINS_URL)
	SETUP_TARGET_FOR_COVERAGE_COBERTURA(testrunner_coverage testrunner "testrunner_coverage_results" "-j;${PROCESSOR_COUNT}")
ELSE()
	SETUP_TARGET_FOR_COVERAGE(testrunner_coverage testrunner "testrunner_coverage_results" "-j;${PROCESSOR_COUNT}")
	SETUP_TARGET_FOR_COVERAGE(ogs-gui_coverage ${EXECUTABLE_OUTPUT_PATH}/ogs-gui "ogs-gui_coverage_results")
ENDIF()
