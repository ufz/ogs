# Script-mode file, run via `cmake -P` from the `web` custom target.
# Expects CTEST_COMMAND, SOURCE_DIR and BINARY_DIR to be passed with -D.
execute_process(
    COMMAND ${CTEST_COMMAND} -N --test-dir ${BINARY_DIR}
    OUTPUT_VARIABLE _ctest_output
    COMMAND_ERROR_IS_FATAL ANY
)
string(REGEX MATCH "Total Tests: ([0-9]+)" _ "${_ctest_output}")
if(NOT _)
    message(FATAL_ERROR "Could not determine the CTest count from:\n${_ctest_output}")
endif()
set(OGS_WEB_CTEST_COUNT ${CMAKE_MATCH_1})

configure_file(${SOURCE_DIR}/web/stats.toml.in
               ${BINARY_DIR}/web/data/stats.toml @ONLY)
