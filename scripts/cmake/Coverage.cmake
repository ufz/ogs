if(NOT (OGS_COVERAGE AND PROJECT_IS_TOP_LEVEL))
    return()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Og --coverage")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fprofile-abs-path")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "--coverage")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "--coverage")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "--coverage")
else()
    message(FATAL_ERROR "OGS_COVERAGE requires clang or gcc compiler!")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    execute_process(
        COMMAND xcrun --find gcov OUTPUT_VARIABLE GCOV_EXECUTABLE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_program(LLVM_COV_EXECUTABLE llvm-cov REQUIRED)
    file(CREATE_LINK ${LLVM_COV_EXECUTABLE} ${CMAKE_BINARY_DIR}/gcov SYMBOLIC)
    set(GCOV_EXECUTABLE "${LLVM_COV_EXECUTABLE} gcov")
else() # Assuming gcc for this example
    find_program(GCOV_EXECUTABLE gcov REQUIRED)
endif()
configure_file(scripts/cmake/gcovr.cfg.in gcovr.cfg @ONLY)

find_program(GCOVR_EXECUTABLE NAMES gcovr)
if(NOT GCOVR_EXECUTABLE)
    list(APPEND OGS_PYTHON_PACKAGES "gcovr==6.0")
    set(GCOVR_EXECUTABLE ${LOCAL_VIRTUALENV_BIN_DIR}/gcovr CACHE PATH "" FORCE)
endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/coverage/html)

add_custom_target(
    process_coverage
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running gcovr to process coverage results"
    COMMAND ${GCOVR_EXECUTABLE} --config gcovr.cfg .
)

if(UNIX)
    add_custom_target(clean_coverage find . -name '*.gcda' -delete)
endif()
