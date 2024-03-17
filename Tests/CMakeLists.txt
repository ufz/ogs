set(CMAKE_FOLDER "Testing")
# ctest dependencies
foreach(tool vtkdiff xdmfdiff)
    if(TARGET ${tool})
        list(APPEND test_dependencies ${tool})
    endif()
endforeach()
add_dependencies(ctest ${test_dependencies})
add_dependencies(ctest-large ${test_dependencies})

# testrunner ###

# VS2012 doesn't support correctly the tuples yet See
# http://code.google.com/p/googletest/issues/detail?id=412
if(MSVC)
    add_definitions(/D_VARIADIC_MAX=10)
endif()

get_source_files(TEST_SOURCES)
append_source_files(TEST_SOURCES BaseLib)
append_source_files(TEST_SOURCES FileIO)
append_source_files(TEST_SOURCES GeoLib)
append_source_files(TEST_SOURCES GeoLib/IO)
append_source_files(TEST_SOURCES MaterialLib)
append_source_files(TEST_SOURCES MaterialLib/SolidModels)
append_source_files(TEST_SOURCES MathLib)
append_source_files(TEST_SOURCES MeshLib)
append_source_files(TEST_SOURCES MeshToolsLib)
append_source_files(TEST_SOURCES MeshGeoToolsLib)
append_source_files(TEST_SOURCES_NUMLIB NumLib)
# Disable Unity build for some files
if("$ENV{HOSTNAME}" MATCHES "frontend.*")
    execute_process(
        COMMAND pwd -P
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE _current_source_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(_current_source_dir ${CMAKE_CURRENT_SOURCE_DIR})
endif()
set(_unity_disabled
    ${_current_source_dir}/NumLib/TestFe.cpp
    ${_current_source_dir}/NumLib/TestGradShapeFunction.cpp
    ${_current_source_dir}/NumLib/TestODEInt.cpp
    ${_current_source_dir}/NumLib/TestCoordinateMapping.cpp
)

set_source_files_properties(
    ${_unity_disabled} PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
# Keep order of source files, otherwise PETSc non-MPI tests fail
set(TEST_SOURCES ${TEST_SOURCES} ${TEST_SOURCES_NUMLIB})
append_source_files(TEST_SOURCES ParameterLib)

if(NOT OGS_BUILD_WHEEL)
    append_source_files(TEST_SOURCES ProcessLib)
    if(OGS_BUILD_PROCESS_TH2M)
        append_source_files(TEST_SOURCES ProcessLib/TH2M)
    endif()
endif()

append_source_files(TEST_SOURCES ProcessLib/Graph)

if(OGS_BUILD_GUI)
    append_source_files(TEST_SOURCES FileIO_Qt)
endif()

if(TARGET ComponentTransport)
    append_source_files(TEST_SOURCES ProcessLib/ComponentTransport)
endif()

if(OGS_BUILD_PROCESS_RichardsMechanics)
    append_source_files(TEST_SOURCES ProcessLib/RichardsMechanics)
endif()

if(OGS_USE_PETSC)
    list(REMOVE_ITEM TEST_SOURCES NumLib/TestSerialLinearSolver.cpp)
endif()

if(OGS_USE_MFRONT)
    add_subdirectory(MaterialLib/MFront)
    append_source_files(TEST_SOURCES MaterialLib/MFront)
endif()

ogs_add_executable(testrunner ${TEST_SOURCES})

target_sources(testrunner
  PRIVATE
    ProcessLib/Graph/TestGet.cpp
    ProcessLib/Graph/TestApply.cpp
)

target_link_libraries(
    testrunner
    PRIVATE $<$<NOT:$<BOOL:${OGS_BUILD_WHEEL}>>:ApplicationsLib>
            ApplicationsFileIO
            autocheck
            gmock
            gtest
            MeshToolsLib
            MeshGeoToolsLib
            MaterialLib
            MathLib
            MeshLib
            MeshToolsLib
            NumLib
            ParameterLib
            $<$<NOT:$<BOOL:${OGS_BUILD_WHEEL}>>:ProcessLib>
            $<$<NOT:$<BOOL:${OGS_BUILD_WHEEL}>>:Processes>
            TestInfoLib
            VTK::FiltersGeneral
            VTK::FiltersSources
            Boost::boost # TODO move OGSMFrontThermodynamicForcesView to
                         # MathLib?
            $<$<TARGET_EXISTS:LIE>:LIE>
            $<$<TARGET_EXISTS:TH2M>:TH2M>
            $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
            $<$<TARGET_EXISTS:SwmmInterface>:SwmmInterface>
            $<$<TARGET_EXISTS:InSituLib>:InSituLib>
            $<$<TARGET_EXISTS:petsc>:petsc>
    PUBLIC $<$<BOOL:${OGS_USE_MFRONT}>:OgsMFrontBehaviourForUnitTests>
)

target_compile_definitions(
    testrunner
    PUBLIC # Enabled elements
           OGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM}
           OGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER}
           $<$<BOOL:${OGS_ENABLE_ELEMENT_SIMPLEX}>:OGS_ENABLE_ELEMENT_SIMPLEX>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_CUBOID}>:OGS_ENABLE_ELEMENT_CUBOID>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PRISM}>:OGS_ENABLE_ELEMENT_PRISM>
           $<$<BOOL:${OGS_ENABLE_ELEMENT_PYRAMID}>:OGS_ENABLE_ELEMENT_PYRAMID>
)

if(OGS_BUILD_GUI)
    target_compile_definitions(testrunner PUBLIC OGS_BUILD_GUI)
    target_link_libraries(
        testrunner
        PRIVATE GitInfoLib
                Qt5::Core
                Qt5::Gui
                Qt5::Xml
                Qt5::Network
                VtkVisFilter
                QtDataView
                VTK::IOImage
    )
    if(GEOTIFF_FOUND)
        target_link_libraries(testrunner PRIVATE ${GEOTIFF_LIBRARIES})
    endif()
endif()

target_precompile_headers(
    testrunner
    PRIVATE
    [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
    <Eigen/Core>
    <gtest/gtest.h>
)

# cmake-format: off
# Add make-target tests which runs the testrunner
if(DEFINED ENV{CI} AND NOT OGS_COVERAGE)
    set(OGS_CI_TESTRUNNER_REPEAT 3 CACHE STRING
        "The number of testrunner repeats for the tests target")
    set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        --gtest_shuffle --gtest_repeat=${OGS_CI_TESTRUNNER_REPEAT})
endif()
set(TESTRUNNER_ADDITIONAL_ARGUMENTS ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
    -l warn
    --gtest_output=xml:./testrunner.xml)

add_custom_target(tests-cleanup ${CMAKE_COMMAND} -E remove -f testrunner.xml)

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
else()
    add_custom_target(tests
        $<TARGET_FILE:testrunner> ${TESTRUNNER_ADDITIONAL_ARGUMENTS}
        DEPENDS testrunner tests-cleanup
    )
endif()
# cmake-format: on

unset(CMAKE_FOLDER)
