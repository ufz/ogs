# Returns the current subdirectory in the sources directory.
macro(GET_CURRENT_SOURCE_SUBDIRECTORY CURRENT_SOURCE_SUBDIRECTORY)
    string(REGEX REPLACE ".*/([^/]*)" "\\1" REGEX_RESULT "${CMAKE_CURRENT_SOURCE_DIR}" )
    set(${CURRENT_SOURCE_SUBDIRECTORY} ${REGEX_RESULT})
endmacro()

# Returns a list of source files (*.h and *.cpp) in SOURCE_FILES and creates a Visual
# Studio folder. A (relative) subdirectory can be passed as second parameter (optional).
macro(GET_SOURCE_FILES SOURCE_FILES)
    if(${ARGC} EQUAL 2)
        set(DIR "${ARGV1}")
    else()
        set(DIR ".")
    endif()

    # Get all files in the directory
    file(GLOB GET_SOURCE_FILES_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.h)
    file(GLOB GET_SOURCE_FILES_TEMPLATES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.tpp)
    file(GLOB GET_SOURCE_FILES_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.cpp)

    set(${SOURCE_FILES} ${GET_SOURCE_FILES_HEADERS} ${GET_SOURCE_FILES_TEMPLATES} ${GET_SOURCE_FILES_SOURCES})
    list(LENGTH ${SOURCE_FILES} NUM_FILES)
    if(${NUM_FILES} EQUAL 0)
        message(FATAL_ERROR "No source files found in ${DIR}")
    endif()

    # Adapt DIR var to backslash syntax of SOURCE_GROUP cmd
    if(${ARGC} EQUAL 2)
        string(REPLACE "/" "\\\\" DIR ${DIR})
        set(DIR "\\${DIR}")
    else()
        set(DIR "")
    endif()

    GET_CURRENT_SOURCE_SUBDIRECTORY(DIRECTORY)
    source_group("${DIRECTORY}${DIR}" FILES
        ${GET_SOURCE_FILES_HEADERS}
        ${GET_SOURCE_FILES_SOURCES}
        ${GET_SOURCE_FILES_TEMPLATES})

endmacro()

# Appends a list of source files (*.h and *.cpp) to SOURCE_FILES and creates a Visual
# Studio folder. A (relative) subdirectory can be passed as second parameter (optional).
macro(APPEND_SOURCE_FILES SOURCE_FILES)
    if(${ARGC} EQUAL 2)
        set(DIR "${ARGV1}")
    else()
        set(DIR ".")
    endif()

    GET_SOURCE_FILES(TMP_SOURCES "${DIR}")
    set(${SOURCE_FILES} ${${SOURCE_FILES}} ${TMP_SOURCES})
endmacro()

# Creates one ctest for each googletest found in source files passed as arguments
# number two onwards. Argument one specifies the testrunner executable.
macro(ADD_GOOGLE_TESTS executable)
    foreach(source ${ARGN})
        file(READ "${source}" contents)
        string(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
        foreach(hit ${found_tests})
            string(REGEX REPLACE ".*\\(([A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+)\\).*" "\\1.\\2" test_name ${hit})
            add_test(${test_name} ${executable}  --gtest_output=xml --gtest_filter=${test_name} ${MI3CTestingDir})
        endforeach()
    endforeach()
endmacro()
