# Supply include directories and compiler flags
get_directory_property(INCLUDE_DIRS INCLUDE_DIRECTORIES)
set(CMAKE_REQUIRED_INCLUDES ${INCLUDE_DIRS})
set(CMAKE_REQUIRED_FLAGS "-std=gnu++14")
set(CMAKE_REQUIRED_QUIET TRUE)

# Checks header for standalone compilation
function(check_header_compilation)
    if(NOT OGS_CHECK_HEADER_COMPILATION)
        return()
    endif()
    string(REPLACE "${PROJECT_SOURCE_DIR}/" "" DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    message(STATUS "Checking header compilation for ${DIRECTORY} ...")
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
    file(GLOB_RECURSE FILES *.h)
    foreach(FILE ${FILES})
        # Ignore *-impl.h files
        if("${FILE}" MATCHES ".*-impl.h")
            continue()
        endif()
        check_cxx_source_compiles(
            "
            #include \"${FILE}\"
            int main() { return 0; }
            "
            COMPILES
        )
        if(NOT COMPILES)
            string(REPLACE "${PROJECT_SOURCE_DIR}/" "" FILE_SHORT ${FILE})
            message(STATUS "  Compilation failed for ${FILE_SHORT}")
        endif()
        unset(COMPILES CACHE)
    endforeach()
endfunction()
