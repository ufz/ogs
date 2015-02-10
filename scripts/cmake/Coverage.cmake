include(CodeCoverage)

set(COVERAGE_EXCLUDES
	'/gpfs0/*'
	'/usr/*'
	'${CMAKE_BINARY_DIR}/*'
	'${CMAKE_SOURCE_DIR}/Tests/*'
	'${CMAKE_SOURCE_DIR}/ThirdParty/*'
)

if(LCOV_PATH AND GENHTML_PATH)
	SETUP_TARGET_FOR_COVERAGE(testrunner_coverage testrunner "testrunner_coverage_results" "-j;${PROCESSOR_COUNT}")
	if(OGS_BUILD_GUI)
		SETUP_TARGET_FOR_COVERAGE(DataExplorer_coverage DataExplorer "DataExplorer_coverage_results")
	endif()
else()
	message(STATUS "No lcov coverage report generated because lcov or genhtml was not found.")
endif()

if(PYTHON_EXECUTABLE)
	SETUP_TARGET_FOR_COVERAGE_COBERTURA(testrunner_coverage_cobertura testrunner "testrunner_coverage_cobertura_results" "-j;${PROCESSOR_COUNT}")
else()
	message(STATUS "No cobertura coverage report generated because Python executable was not found.")
endif()
