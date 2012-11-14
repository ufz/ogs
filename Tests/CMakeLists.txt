# VS2012 doesn't support correctly the tuples yet
# See http://code.google.com/p/googletest/issues/detail?id=412
IF( MSVC )
	ADD_DEFINITIONS( /D_VARIADIC_MAX=10 )
ENDIF()

APPEND_SOURCE_FILES(TEST_SOURCES)
APPEND_SOURCE_FILES(TEST_SOURCES BaseLib)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib)
APPEND_SOURCE_FILES(TEST_SOURCES MathLib)

INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Tests/gtest
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
)

ADD_EXECUTABLE (testrunner testrunner.cpp ${TEST_SOURCES})
SET_TARGET_PROPERTIES(testrunner PROPERTIES FOLDER Testing)

TARGET_LINK_LIBRARIES(testrunner
	GTest
	BaseLib
	GeoLib
	MathLib
)

IF(USE_EXTERNAL_LIS)
    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
    TARGET_LINK_LIBRARIES( testrunner ${LIS_LIBRARIES} )
ELSE()
    IF (USE_LIS)
        INCLUDE_DIRECTORIES( ../ThirdParty/lis/include )
        TARGET_LINK_LIBRARIES( testrunner lis )
    ENDIF()
ENDIF ()

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
