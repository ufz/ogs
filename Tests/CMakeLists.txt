# VS2012 doesn't support correctly the tuples yet
# See http://code.google.com/p/googletest/issues/detail?id=412
IF( MSVC )
	ADD_DEFINITIONS( /D_VARIADIC_MAX=10 )
ENDIF()

APPEND_SOURCE_FILES(TEST_SOURCES)
APPEND_SOURCE_FILES(TEST_SOURCES BaseLib)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib)
APPEND_SOURCE_FILES(TEST_SOURCES MathLib)
APPEND_SOURCE_FILES(TEST_SOURCES MeshLib)

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

IF (LIS_FOUND)
	INCLUDE_DIRECTORIES(${LIS_INCLUDE_DIR})
ENDIF()

ADD_EXECUTABLE (testrunner testrunner.cpp ${TEST_SOURCES})
SET_TARGET_PROPERTIES(testrunner PROPERTIES FOLDER Testing)

TARGET_LINK_LIBRARIES(testrunner
	GTest
	BaseLib
	GeoLib
	MathLib
	MeshLib
	logog
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
)

# Add make-target test which runs the testrunner
# This should override CTest's predefined test-target but it does not
ADD_CUSTOM_TARGET(test
	$<TARGET_FILE:testrunner>
	DEPENDS testrunner
)

# Testing
#ENABLE_TESTING()

# Creates one ctest entry for every googletest
#ADD_GOOGLE_TESTS ( ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/testrunner ${TEST_SOURCES})
