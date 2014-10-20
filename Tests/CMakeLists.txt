# VS2012 doesn't support correctly the tuples yet
# See http://code.google.com/p/googletest/issues/detail?id=412
IF( MSVC )
	ADD_DEFINITIONS( /D_VARIADIC_MAX=10 )
ENDIF()

APPEND_SOURCE_FILES(TEST_SOURCES)
APPEND_SOURCE_FILES(TEST_SOURCES AssemblerLib)
APPEND_SOURCE_FILES(TEST_SOURCES BaseLib)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib)
APPEND_SOURCE_FILES(TEST_SOURCES MathLib)
APPEND_SOURCE_FILES(TEST_SOURCES MeshLib)
APPEND_SOURCE_FILES(TEST_SOURCES NumLib)

IF (QT4_FOUND)
	APPEND_SOURCE_FILES(TEST_SOURCES FileIO)
ENDIF()

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/AssemblerLib
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/MeshGeoToolsLib
	${CMAKE_SOURCE_DIR}/NumLib
	${CMAKE_BINARY_DIR}/BaseLib
)

ADD_EXECUTABLE (testrunner testrunner.cpp ${TEST_SOURCES})
SET_TARGET_PROPERTIES(testrunner PROPERTIES FOLDER Testing)
IF(TARGET Eigen)
	ADD_DEPENDENCIES(testrunner Eigen)
ENDIF()

TARGET_LINK_LIBRARIES(testrunner
	GTest
	ApplicationsLib
	AssemblerLib
	BaseLib
	FileIO
	GeoLib
	MathLib
	MeshLib
	MeshGeoToolsLib
	NumLib
	logog
	${BOOST_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
)
IF (OGS_USE_PETSC)
	TARGET_LINK_LIBRARIES( testrunner  ${PETSC_LIBRARIES})
ENDIF (OGS_USE_PETSC)

IF (OGS_USE_MPI)
	TARGET_LINK_LIBRARIES( testrunner  ${MPI_CXX_LIBRARIES})
ENDIF ()

IF(OGS_BUILD_GUI)
	TARGET_LINK_LIBRARIES(testrunner
		QtDataView
		VtkVis
		${VTK_LIBRARIES}
	)
ENDIF()

IF (QT4_FOUND)
	TARGET_LINK_LIBRARIES(testrunner ${QT_LIBRARIES})
ENDIF()


# Add make-target test which runs the testrunner
# This should override CTest's predefined test-target but it does not
IF (OGS_USE_PETSC)
	ADD_CUSTOM_TARGET(tests
		mpirun -np 1 $<TARGET_FILE:testrunner> --gtest_filter=-MPITest*
		DEPENDS testrunner
	)
	ADD_CUSTOM_TARGET(tests_mpi
		mpirun -np 3 $<TARGET_FILE:testrunner>  --gtest_filter=MPITest*
		DEPENDS testrunner
	)
ELSE ()
	ADD_CUSTOM_TARGET(tests
		$<TARGET_FILE:testrunner>
		DEPENDS testrunner
	)
ENDIF ()

# Testing
#ENABLE_TESTING()

# Creates one ctest entry for every googletest
#ADD_GOOGLE_TESTS ( ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/testrunner ${TEST_SOURCES})
