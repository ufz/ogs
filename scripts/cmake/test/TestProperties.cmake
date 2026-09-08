# Shared test registration properties.

if(DEFINED OGS_TEST_PROPERTIES_INCLUDED)
    return()
endif()
set(OGS_TEST_PROPERTIES_INCLUDED TRUE)

set(OGS_CTEST_OMP_THREADS 4
    CACHE STRING "Number of threads used by OpenMP ctest variants."
)

# Resolve the runtime policy shared by AddTest, OgsTest, and NotebookTest.
# Outputs are the resolved runtime, an optional timeout, the large-test flag,
# and the skip flag, in that order.

# cmake-lint: disable=R0913
function(
    ogs_resolve_test_runtime
    runtime
    large_runtime
    max_runtime
    out_runtime
    out_timeout
    out_is_large
    out_skip
)
    if("${runtime}" STREQUAL "")
        set(runtime_value 1)
    else()
        set(runtime_value "${runtime}")
    endif()

    set(timeout_value)
    if(runtime_value GREATER 750)
        math(EXPR timeout_value "${runtime_value} * 2")
    endif()

    set(skip_value FALSE)
    if(NOT "${max_runtime}" STREQUAL "")
        if(runtime_value GREATER ${max_runtime})
            set(skip_value TRUE)
        endif()
    endif()

    if(runtime_value GREATER ${large_runtime})
        set(is_large_value TRUE)
    else()
        set(is_large_value FALSE)
    endif()

    set(${out_runtime} "${runtime_value}" PARENT_SCOPE)
    set(${out_timeout} "${timeout_value}" PARENT_SCOPE)
    set(${out_is_large} "${is_large_value}" PARENT_SCOPE)
    set(${out_skip} "${skip_value}" PARENT_SCOPE)
endfunction()

# Add OpenMP properties while retaining the test's existing environment.
function(ogs_set_omp_test_properties test_name base_processors labels)
    get_test_property(${test_name} ENVIRONMENT existing_test_environment)
    if(NOT existing_test_environment)
        set(existing_test_environment "")
    endif()
    math(EXPR overall_processors
         "${OGS_CTEST_OMP_THREADS} * ${base_processors}"
    )
    set_tests_properties(
        ${test_name}
        PROPERTIES
            ENVIRONMENT
            "OGS_ASM_THREADS=${OGS_CTEST_OMP_THREADS};${existing_test_environment}"
            PROCESSORS
            ${overall_processors}
            LABELS
            "${labels};omp"
    )
endfunction()
