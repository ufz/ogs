if(NOT CPPCHECK_TOOL_PATH)
    return()
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CPPCHECK_PARALLEL 4)
if(DEFINED ENV{CMAKE_BUILD_PARALLEL_LEVEL})
    set(CPPCHECK_PARALLEL $ENV{CMAKE_BUILD_PARALLEL_LEVEL})
elseif(DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
    set(CPPCHECK_PARALLEL ${CMAKE_BUILD_PARALLEL_LEVEL})
endif()
if(DEFINED CPM_SOURCE_CACHE)
    set(_cpm_dir ${CPM_SOURCE_CACHE})
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
    set(_cpm_dir $ENV{CPM_SOURCE_CACHE})
else()
    set(_cpm_dir ${PROJECT_BINARY_DIR}/_deps)
endif()
set(_last_sed_expression [['x;${s/,$//;p;x;};1d']])

configure_file(
    ${PROJECT_SOURCE_DIR}/scripts/test/cppcheck.in.sh
    ${PROJECT_BINARY_DIR}/cppcheck.sh
)

if(DEFINED ENV{NUM_THREADS})
    set(CPPCHECK_THREADS -j $ENV{NUM_THREADS})
endif()

add_custom_target(
    cppcheck
    COMMAND ${BASH_TOOL_PATH} cppcheck.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)
