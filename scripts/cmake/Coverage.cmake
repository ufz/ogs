INCLUDE(CodeCoverage)

SET(COVERAGE_EXCLUDES
	'/usr/*'
	'${CMAKE_BINARY_DIR}/*'
	'${CMAKE_SOURCE_DIR}/tests/*'
	'${CMAKE_SOURCE_DIR}/BaseLib/zlib/*'
	'${CMAKE_SOURCE_DIR}/BaseLib/logog/*'
	'${CMAKE_SOURCE_DIR}/BaseLib/RapidXML/*'
)

IF(JENKINS_URL)
	SETUP_TARGET_FOR_COVERAGE_COBERTURA(ctest_coverage ctest "ctest_coverage_results" "-j;${PROCESSOR_COUNT}")
ELSE()
	SETUP_TARGET_FOR_COVERAGE(ctest_coverage ctest "ctest_coverage_results" "-j;${PROCESSOR_COUNT}")
	SETUP_TARGET_FOR_COVERAGE(ogs-gui_coverage ${EXECUTABLE_OUTPUT_PATH}/ogs-gui "ogs-gui_coverage_results")
ENDIF()