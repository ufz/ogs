# VS2012 doesn't support correctly the tuples yet
# See http://code.google.com/p/googletest/issues/detail?id=412
if(MSVC)
	add_definitions(/D_VARIADIC_MAX=10)
endif()

APPEND_SOURCE_FILES(TEST_SOURCES)
APPEND_SOURCE_FILES(TEST_SOURCES AssemblerLib)
APPEND_SOURCE_FILES(TEST_SOURCES BaseLib)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib)
if(NOT OGS_DONT_USE_VTK)
	APPEND_SOURCE_FILES(TEST_SOURCES InSituLib)
endif()
APPEND_SOURCE_FILES(TEST_SOURCES MathLib)
APPEND_SOURCE_FILES(TEST_SOURCES MeshLib)
APPEND_SOURCE_FILES(TEST_SOURCES NumLib)

set(TEST_SOURCES ${TEST_SOURCES} FileIO/TestGLIReader.cpp)
if(QT4_FOUND)
	set(TEST_SOURCES ${TEST_SOURCES} FileIO/TestXmlGmlReader.cpp)
endif()

include_directories(
	${CMAKE_SOURCE_DIR}/AssemblerLib
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/InSituLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/MeshGeoToolsLib
	${CMAKE_SOURCE_DIR}/NumLib
	${CMAKE_BINARY_DIR}/BaseLib
)

add_executable(testrunner testrunner.cpp ${TEST_SOURCES})
set_target_properties(testrunner PROPERTIES FOLDER Testing)
if(TARGET Eigen)
	add_dependencies(testrunner Eigen)
endif()

target_link_libraries(testrunner
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
if(NOT OGS_DONT_USE_VTK)
	target_link_libraries(testrunner
		InSituLib
		${CATALYST_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(testrunner)
endif()
if(OGS_USE_PETSC)
	target_link_libraries(testrunner ${PETSC_LIBRARIES})
endif()

if(OGS_USE_MPI)
	target_link_libraries(testrunner ${MPI_CXX_LIBRARIES})
endif()

if(OGS_BUILD_GUI)
	target_link_libraries(testrunner
		QtDataView
		VtkVis
	)
endif()

if(QT4_FOUND)
	target_link_libraries(testrunner ${QT_LIBRARIES})
endif()

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(testrunner)

# Add make-target test which runs the testrunner
# This should override CTest's predefined test-target but it does not
if(DEFINED ENV{CI})
	set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS} --gtest_shuffle --gtest_repeat=3)
endif()
if(OGS_USE_PETSC)
	add_custom_target(tests
		mpirun -np 1 $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS} --gtest_filter=-MPITest*:*AssemblerLib*
		DEPENDS testrunner
	)
	add_custom_target(tests_mpi
		mpirun -np 3 $<TARGET_FILE:testrunner> --gtest_filter=MPITest*
		DEPENDS testrunner
	)
else()
	add_custom_target(tests
		$<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
		DEPENDS testrunner
	)
endif()

# Creates one ctest entry for every googletest
#ADD_GOOGLE_TESTS ( ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/testrunner ${TEST_SOURCES})
