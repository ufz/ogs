### xdmfdiff ###
if(OGS_USE_XDMF)
    add_subdirectory(xdmfdiff)
endif()

### testrunner ###
include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# VS2012 doesn't support correctly the tuples yet See
# http://code.google.com/p/googletest/issues/detail?id=412
if(MSVC)
    add_definitions(/D_VARIADIC_MAX=10)
endif()

append_source_files(TEST_SOURCES)
append_source_files(TEST_SOURCES BaseLib)
append_source_files(TEST_SOURCES FileIO)
append_source_files(TEST_SOURCES GeoLib)
append_source_files(TEST_SOURCES GeoLib/IO)
append_source_files(TEST_SOURCES MaterialLib)
append_source_files(TEST_SOURCES MathLib)
append_source_files(TEST_SOURCES MeshLib)
append_source_files(TEST_SOURCES MeshGeoToolsLib)
append_source_files(TEST_SOURCES_NUMLIB NumLib)
# Disable Unity build for NumLib tests
set_source_files_properties(${TEST_SOURCES_NUMLIB} PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
# Keep order of source files, otherwise PETSc non-MPI tests fail
set(TEST_SOURCES ${TEST_SOURCES} ${TEST_SOURCES_NUMLIB})
append_source_files(TEST_SOURCES ParameterLib)
append_source_files(TEST_SOURCES ProcessLib)

if(OGS_BUILD_GUI)
    append_source_files(TEST_SOURCES FileIO_Qt)
endif()

if(OGS_BUILD_SWMM)
    append_source_files(TEST_SOURCES FileIO_SWMM)
endif()

if(OGS_BUILD_PROCESS_LIE)
    append_source_files(TEST_SOURCES ProcessLib/LIE)
endif()

if(OGS_USE_PETSC)
    list(REMOVE_ITEM TEST_SOURCES NumLib/TestSerialLinearSolver.cpp)
endif()

add_executable(testrunner ${TEST_SOURCES})
set_target_properties(testrunner PROPERTIES FOLDER Testing)

target_link_libraries(testrunner
    PRIVATE
        ApplicationsFileIO
        gtest
        MeshGeoToolsLib
        MaterialLib
        MeshLib
        NumLib
        ParameterLib
        ProcessLib
        TestInfoLib
        Threads::Threads
        ${VTK_LIBRARIES}
        $<$<TARGET_EXISTS:LIE>:LIE>
        $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
        $<$<TARGET_EXISTS:SwmmInterface>:SwmmInterface>
        $<$<TARGET_EXISTS:InSituLib>:InSituLib>
        $<$<BOOL:${OGS_USE_PETSC}>:${PETSC_LIBRARIES}>
)

if(OGS_BUILD_GUI)
    target_link_libraries(testrunner
        PRIVATE
            GitInfoLib
            Qt5::Core
            Qt5::Gui
            Qt5::Xml
            Qt5::Network
            VtkVis
            QtDataView
    )
    if(GEOTIFF_FOUND)
        target_link_libraries(testrunner PRIVATE ${GEOTIFF_LIBRARIES})
    endif()
endif()

# cmake-format: off
# Add make-target tests which runs the testrunner
if(IS_CI AND NOT OGS_COVERAGE)
    set(OGS_CI_TESTRUNNER_REPEAT 3 CACHE STRING
        "The number of testrunner repeats for the tests target")
    set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        --gtest_shuffle --gtest_repeat=${OGS_CI_TESTRUNNER_REPEAT})
endif()
set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
    -l warn
    --gtest_output=xml:./testrunner.xml)

add_custom_target(tests-cleanup ${CMAKE_COMMAND} -E remove -f testrunner.xml)
set_target_properties(tests-cleanup PROPERTIES FOLDER Testing)

if(OGS_USE_PETSC)
    if("${HOSTNAME}" MATCHES "frontend.*")
        set(MPIRUN_ARGS --mca btl_openib_allow_ib 1)
    endif()
    set(TEST_FILTER_MPI --gtest_filter=-MPITest*)
    add_custom_target(tests
        mpirun ${MPIRUN_ARGS} -np 1 $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS} ${TEST_FILTER_MPI}
        DEPENDS testrunner tests-cleanup
    )
    add_custom_target(tests_mpi
        mpirun ${MPIRUN_ARGS} -np 3 $<TARGET_FILE:testrunner> --gtest_filter=MPITest*
        DEPENDS testrunner
    )
    set_target_properties(tests_mpi PROPERTIES FOLDER Testing)
else()
    add_custom_target(tests
        $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        DEPENDS testrunner tests-cleanup
    )
endif()
# cmake-format: on

set_target_properties(tests PROPERTIES FOLDER Testing)

# Creates one ctest entry for every googletest
# ~~~
# ADD_GOOGLE_TESTS (${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/testrunner ${TEST_SOURCES})
# ~~~
