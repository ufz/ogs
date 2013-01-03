# Returns the current subdirectory in the sources directory.
MACRO(GET_CURRENT_SOURCE_SUBDIRECTORY CURRENT_SOURCE_SUBDIRECTORY)
	STRING(REGEX REPLACE ".*/([^/]*)" "\\1" REGEX_RESULT "${CMAKE_CURRENT_SOURCE_DIR}" )
	SET(${CURRENT_SOURCE_SUBDIRECTORY} ${REGEX_RESULT})
ENDMACRO()

# Returns a list of source files (*.h and *.cpp) in SOURCE_FILES and creates a Visual
# Studio folder. A (relative) subdirectory can be passed as second parameter (optional).
MACRO(GET_SOURCE_FILES SOURCE_FILES)

	IF(${ARGC} EQUAL 2)
		SET(DIR "${ARGV1}")
	ELSE()
		SET(DIR ".")
	ENDIF()

	# Get all files in the directory
	FILE(GLOB GET_SOURCE_FILES_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.h)
	FILE(GLOB GET_SOURCE_FILES_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.tpp)
	FILE(GLOB GET_SOURCE_FILES_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.cpp)

	SET(${SOURCE_FILES} ${GET_SOURCE_FILES_HEADERS} ${GET_SOURCE_FILES_SOURCES})

	# Adapt DIR var to backslash syntax of SOURCE_GROUP cmd
	IF(${ARGC} EQUAL 2)
		STRING(REPLACE "/" "\\\\" DIR ${DIR})
		SET(DIR "\\${DIR}")
	ELSE()
		SET(DIR "")
	ENDIF()

	GET_CURRENT_SOURCE_SUBDIRECTORY(DIRECTORY)
	SOURCE_GROUP( "${DIRECTORY}${DIR}" FILES
		${GET_SOURCE_FILES_HEADERS}
		${GET_SOURCE_FILES_SOURCES})

ENDMACRO()

# Appends a list of source files (*.h and *.cpp) to SOURCE_FILES and creates a Visual
# Studio folder. A (relative) subdirectory can be passed as second parameter (optional).
MACRO(APPEND_SOURCE_FILES SOURCE_FILES)
	GET_SOURCE_FILES(TMP_SOURCES "${ARGV}")
	SET(${SOURCE_FILES} ${${SOURCE_FILES}} ${TMP_SOURCES})
ENDMACRO()

# Creates one ctest for each googletest found in source files passed as arguments
# number two onwards. Argument one specifies the testrunner executable.
MACRO(ADD_GOOGLE_TESTS executable)
	FOREACH ( source ${ARGN} )
		FILE(READ "${source}" contents)
		STRING(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
		FOREACH(hit ${found_tests})
			STRING(REGEX REPLACE ".*\\(([A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+)\\).*" "\\1.\\2" test_name ${hit})
			ADD_TEST(${test_name} ${executable}  --gtest_output=xml --gtest_filter=${test_name} ${MI3CTestingDir})
			# message ("Adding test: ${test_name}")
		ENDFOREACH(hit)
	ENDFOREACH()
ENDMACRO()
