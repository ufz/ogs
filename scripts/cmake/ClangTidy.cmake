if(NOT PROJECT_IS_TOP_LEVEL)
    return()
endif()

find_program(CLANG_TIDY_EXECUTABLE clang-tidy)
find_program(RUN_CLANG_TIDY_EXECUTABLE run-clang-tidy)
if(CLANG_TIDY_EXECUTABLE AND RUN_CLANG_TIDY_EXECUTABLE)
    add_custom_target(
        run-clang-tidy
        COMMAND ${RUN_CLANG_TIDY_EXECUTABLE} -clang-tidy-binary
                ${CLANG_TIDY_EXECUTABLE} -p ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        USES_TERMINAL
    )
    # Disable clang-tidy in cpm build dirs, hack: Disable all checks, but one
    # check has to be enabled.
    file(WRITE ${PROJECT_BINARY_DIR}/_deps/.clang-tidy
         "Checks: '-*,boost-use-to-string'"
    )
endif()
