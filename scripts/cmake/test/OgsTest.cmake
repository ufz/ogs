function (OgsTest)
    if(NOT BUILD_TESTING OR NOT OGS_BUILD_CLI)
        return()
    endif()
    set(options LARGE)
    set(oneValueArgs PROJECTFILE RUNTIME)
    set(multiValueArgs WRAPPER)
    cmake_parse_arguments(OgsTest "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    get_filename_component(OgsTest_DIR "${OgsTest_PROJECTFILE}" DIRECTORY)
    get_filename_component(OgsTest_NAME "${OgsTest_PROJECTFILE}" NAME)
    get_filename_component(OgsTest_NAME_WE "${OgsTest_PROJECTFILE}" NAME_WE)

    if (OgsTest_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unparsed argument(s) '${OgsTest_UNPARSED_ARGUMENTS}' to OgsTest call.")
    endif()

    if (NOT DEFINED OgsTest_RUNTIME)
        set(OgsTest_RUNTIME 1)
    endif()

    set(OgsTest_SOURCE_DIR "${Data_SOURCE_DIR}/${OgsTest_DIR}")
    set(OgsTest_BINARY_DIR "${Data_BINARY_DIR}/${OgsTest_DIR}")
    file(MAKE_DIRECTORY ${OgsTest_BINARY_DIR})
    file(TO_NATIVE_PATH "${OgsTest_BINARY_DIR}" OgsTest_BINARY_DIR_NATIVE)

    set(TEST_NAME "ogs-${OgsTest_DIR}/${OgsTest_NAME_WE}")
    # Add wrapper postfix (-mpi for mpirun).
    if (OgsTest_WRAPPER)
        string(REGEX MATCH "^[^ ]+" WRAPPER ${OgsTest_WRAPPER})
        if (WRAPPER STREQUAL "mpirun")
            set(TEST_NAME "${TEST_NAME}-mpi")
        endif()
    endif()
    set(BATS_FILENAME benchmarks.bats)
    # Add -LARGE tag.
    if (${OgsTest_LARGE})
        set(TEST_NAME "${TEST_NAME}-LARGE")
        set(BATS_FILENAME benchmarks-large.bats)
    endif()

    add_test(
        NAME ${TEST_NAME}
        WORKING_DIRECTORY "${OgsTest_BINARY_DIR}"
        COMMAND ${OgsTest_WRAPPER} $<TARGET_FILE:ogs> -r ${OgsTest_SOURCE_DIR} ${OgsTest_SOURCE_DIR}/${OgsTest_NAME})
    # For debugging:
    #message("Adding test with
    #    NAME ${TEST_NAME}
    #    WORKING_DIRECTORY ${OgsTest_BINARY_DIR}
    #    COMMAND ${OgsTest_WRAPPER} $<TARGET_FILE:ogs> -r ${OgsTest_SOURCE_DIR} ${OgsTest_SOURCE_DIR}/${OgsTest_NAME})

    # bats container benchmark runner
    file(APPEND ${PROJECT_BINARY_DIR}/${BATS_FILENAME} "\
@test \"benchmark - ${TEST_NAME}\" {\n\
  rm -r $OUT/${OgsTest_DIR} | true\n\
  mkdir -p $OUT/${OgsTest_DIR}\n\
  run singularity exec $SIF scif run ogs -o $OUT/${OgsTest_DIR} -r $SRC/${OgsTest_DIR} $SRC/${OgsTest_DIR}/${OgsTest_NAME}\n\
  [ \"$status\" -eq 0 ]\n\
}\n\n\
")

    set_tests_properties(${TEST_NAME} PROPERTIES
        ENVIRONMENT VTKDIFF_EXE=$<TARGET_FILE:vtkdiff>
        COST ${OgsTest_RUNTIME})

    if(TARGET ${OgsTest_EXECUTABLE})
        add_dependencies(ctest ${OgsTest_EXECUTABLE})
        add_dependencies(ctest-large ${OgsTest_EXECUTABLE})
    endif()
endfunction()
