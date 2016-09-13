include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# VS2012 doesn't support correctly the tuples yet
# See http://code.google.com/p/googletest/issues/detail?id=412
if(MSVC)
    add_definitions(/D_VARIADIC_MAX=10)
endif()

APPEND_SOURCE_FILES(TEST_SOURCES)
APPEND_SOURCE_FILES(TEST_SOURCES ApplicationsLib)
APPEND_SOURCE_FILES(TEST_SOURCES BaseLib)
APPEND_SOURCE_FILES(TEST_SOURCES FileIO)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib)
APPEND_SOURCE_FILES(TEST_SOURCES GeoLib/IO)
APPEND_SOURCE_FILES(TEST_SOURCES MaterialLib)
APPEND_SOURCE_FILES(TEST_SOURCES MathLib)
APPEND_SOURCE_FILES(TEST_SOURCES MeshLib)
APPEND_SOURCE_FILES(TEST_SOURCES MeshGeoToolsLib)
APPEND_SOURCE_FILES(TEST_SOURCES NumLib)
APPEND_SOURCE_FILES(TEST_SOURCES MaterialLib)

if(QT4_FOUND)
    APPEND_SOURCE_FILES(TEST_SOURCES FileIO_Qt)
endif()

if(OGS_USE_PETSC OR OGS_USE_MPI)
    list(REMOVE_ITEM TEST_SOURCES NumLib/TestSerialLinearSolver.cpp)
endif()

add_executable(testrunner ${TEST_SOURCES})
set_target_properties(testrunner PROPERTIES FOLDER Testing)

target_link_libraries(testrunner
    ApplicationsLib
    ApplicationsFileIO
    GTest
    MeshGeoToolsLib
    MaterialLib
    MeshLib
    NumLib
    Threads::Threads
)
ADD_VTK_DEPENDENCY(testrunner)

if(OGS_USE_PETSC)
    target_link_libraries(testrunner ${PETSC_LIBRARIES})
endif()

if(OGS_USE_MPI)
    target_link_libraries(testrunner ${MPI_CXX_LIBRARIES})
endif()

if(OGS_BUILD_GUI)
    target_link_libraries(testrunner
        QtDataView
        QtStratView
        VtkVis
    )
endif()

if(QT4_FOUND)
    target_link_libraries(testrunner Qt4::QtCore Qt4::QtGui Qt4::QtXml Qt4::QtNetwork)
    if(CMAKE_CROSSCOMPILING)
        target_link_libraries(testrunner
            ${QT_XML_DEPS_LIBRARIES}
            ${QT_GUI_DEPS_LIBRARIES}
            ${QT_NETWORK_DEPS_LIBRARIES})
    endif()
endif()

ADD_VTK_DEPENDENCY(testrunner)

# Add make-target tests which runs the testrunner
if(DEFINED ENV{CI})
    set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        --gtest_shuffle --gtest_repeat=3)
endif()
set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
    -l warn
    --gtest_output=xml:./testrunner.xml)

add_custom_target(tests-cleanup ${CMAKE_COMMAND} -E remove testrunner.xml)

if(OGS_USE_PETSC)
    set(TEST_FILTER_MPI --gtest_filter=-MPITest_Math.*)
    add_custom_target(tests
        mpirun -np 1 $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS} ${TEST_FILTER_MPI}
        DEPENDS testrunner tests-cleanup
    )
    add_custom_target(tests_mpi
        mpirun -np 3 $<TARGET_FILE:testrunner> --gtest_filter=MPITest*
        DEPENDS testrunner
    )
else()
    add_custom_target(tests
        $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        DEPENDS testrunner tests-cleanup
    )
endif()

# Creates one ctest entry for every googletest
#ADD_GOOGLE_TESTS ( ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/testrunner ${TEST_SOURCES})
