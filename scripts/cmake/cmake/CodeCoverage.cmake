# - Enable Code Coverage
#
# 2012-01-31, Lars Bilke
#
# USAGE:
# 1. Copy this file into your cmake modules path
# 2. Add the following line to your CMakeLists.txt:
#      include(CodeCoverage)
#
# 3. set(COVERAGE_EXCLUDES 'dir1/*' 'dir2/*')
# 4. Use the function SETUP_TARGET_FOR_COVERAGE to create a custom make target
#    which runs your test executable and produces a lcov code coverage report.
#

# Check prereqs
find_program( GCOV_PATH gcov )
find_program( LCOV_PATH lcov )
find_program( GENHTML_PATH genhtml )
find_program( GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/scripts/test)

if(NOT GCOV_PATH)
	message(FATAL_ERROR "gcov not found! Aborting...")
endif() # NOT GCOV_PATH

if(NOT CMAKE_COMPILER_IS_GNUCXX AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message(FATAL_ERROR "Compiler is not GNU gcc! Aborting...")
endif() # NOT CMAKE_COMPILER_IS_GNUCXX

if ( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
  message( WARNING "Code coverage results with an optimised (non-Debug) build may be misleading" )
endif() # NOT CMAKE_BUILD_TYPE STREQUAL "Debug"


# Setup compiler options
add_definitions(-fprofile-arcs -ftest-coverage)

if(CMAKE_COMPILER_IS_GNUCXX)
	link_libraries(gcov)
else()
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# Param _targetname     The name of new the custom make target
# Param _testrunner     The name of the target which runs the tests
# Param _outputname     lcov output is generated as _outputname.info
#                       HTML report is generated in _outputname/index.html
# Optional fourth parameter is passed as arguments to _testrunner
#   Pass them in list form, e.g.: "-j;2" for -j 2
function(SETUP_TARGET_FOR_COVERAGE _targetname _testrunner _outputname)

	if(NOT LCOV_PATH)
		message(FATAL_ERROR "lcov not found! Aborting...")
	endif() # NOT LCOV_PATH

	if(NOT GENHTML_PATH)
		message(FATAL_ERROR "genhtml not found! Aborting...")
	endif() # NOT GENHTML_PATH

	# Setup target
	add_custom_target(${_targetname}

		# Cleanup lcov
		${LCOV_PATH} --directory . --zerocounters

		# Run tests
		COMMAND ${_testrunner} ${ARGV3}

		# Capturing lcov counters and generating report
		COMMAND ${LCOV_PATH} --directory . --capture --output-file ${_outputname}.info
		COMMAND ${LCOV_PATH} --remove ${_outputname}.info ${COVERAGE_EXCLUDES} --output-file ${_outputname}.info.cleaned
		COMMAND ${GENHTML_PATH} -o ${_outputname} ${_outputname}.info.cleaned
		COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned

		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		DEPENDS ${_testrunner}
		COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
	)

	# Show info where to find the report
	add_custom_command(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Open ./${_outputname}/index.html in your browser to view the coverage report."
	)

endfunction() # SETUP_TARGET_FOR_COVERAGE

# Param _targetname     The name of new the custom make target
# Param _testrunner     The name of the target which runs the tests
# Param _outputname     cobertura output is generated as _outputname.xml
# Optional fourth parameter is passed as arguments to _testrunner
#   Pass them in list form, e.g.: "-j;2" for -j 2
function(SETUP_TARGET_FOR_COVERAGE_COBERTURA _targetname _testrunner _outputname)

	if(NOT PYTHON_EXECUTABLE)
		message(FATAL_ERROR "Python not found! Aborting...")
	endif() # NOT PYTHON_EXECUTABLE

	if(NOT GCOVR_PATH)
		message(FATAL_ERROR "gcovr not found! Aborting...")
	endif() # NOT GCOVR_PATH

	# Combine excludes to several -e arguments
	set(COBERTURA_EXCLUDES "")
	foreach(EXCLUDE ${COVERAGE_EXCLUDES})
		set(COBERTURA_EXCLUDES "-e ${EXCLUDE} ${COBERTURA_EXCLUDES}")
	endforeach()

	if(DEFINED ENV{JENKINS_URL})
		# Get relative paths
		string(REPLACE "$ENV{WORKSPACE}/" "" CMAKE_RELATIVE_SOURCE_DIR ${CMAKE_SOURCE_DIR})
		string(REPLACE "/gpfs0" "" CMAKE_RELATIVE_SOURCE_DIR ${CMAKE_RELATIVE_SOURCE_DIR})
		string(REPLACE "/gpfs1" "" CMAKE_RELATIVE_SOURCE_DIR ${CMAKE_RELATIVE_SOURCE_DIR})

		add_custom_target(${_targetname}

			# Run tests
			${_testrunner} ${ARGV3}

			# Running gcovr
			COMMAND ${GCOVR_PATH} -x -d -r ${CMAKE_RELATIVE_SOURCE_DIR} ${COBERTURA_EXCLUDES}
				-o ${CMAKE_BINARY_DIR}/${_outputname}.xml
			WORKING_DIRECTORY $ENV{WORKSPACE}
			DEPENDS ${_testrunner}
			COMMENT "Running gcovr to produce Cobertura code coverage report for Jenkins."
		)
	else()
		add_custom_target(${_targetname}

			# Run tests
			${_testrunner} ${ARGV3}

			# Running gcovr
			COMMAND ${GCOVR_PATH} -x -r ${CMAKE_SOURCE_DIR} ${COBERTURA_EXCLUDES}
				-o ${_outputname}.xml
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			DEPENDS ${_testrunner}
			COMMENT "Running gcovr to produce Cobertura code coverage report."
		)
	endif()

	# Show info where to find the report
	add_custom_command(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Cobertura code coverage report saved in ${_outputname}.xml."
	)

endfunction() # SETUP_TARGET_FOR_COVERAGE_COBERTURA
