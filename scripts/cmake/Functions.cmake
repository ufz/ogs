# Returns the current subdirectory in the sources directory.
MACRO(GET_CURRENT_SOURCE_SUBDIRECTORY CURRENT_SOURCE_SUBDIRECTORY)
	STRING(REGEX REPLACE ".*/([^/]*)" "\\1" REGEX_RESULT "${CMAKE_CURRENT_SOURCE_DIR}" )
	SET(${CURRENT_SOURCE_SUBDIRECTORY} ${REGEX_RESULT})
ENDMACRO()

# Returns a list of source files (*.h and *.cpp) in SOURCE_FILES and creates a Visual
# Studio folder. A (relative) directory can be passed as second parameter (optional).
MACRO(GET_SOURCE_FILES SOURCE_FILES)
	
	IF(ARGC EQUAL 2)
		SET(DIR "${arg2}/")
	ELSE()
		SET(DIR "")
	ENDIF()
	
	# Get all files in the directory
	FILE(GLOB GET_SOURCE_FILES_HEADERS ${DIR}*.h)
	FILE(GLOB GET_SOURCE_FILES_SOURCES ${DIR}*.cpp)

	SET(${SOURCE_FILES} ${GET_SOURCE_FILES_HEADERS} ${GET_SOURCE_FILES_SOURCES})
	
	GET_CURRENT_SOURCE_SUBDIRECTORY(DIRECTORY)
	MESSAGE("Curr: ${DIRECTORY}")
	SOURCE_GROUP( ${DIRECTORY} FILES ${SOURCE_FILES})
	
ENDMACRO()