# Supply include directories and compiler flags
get_directory_property(INCLUDE_DIRS INCLUDE_DIRECTORIES)
set(CMAKE_REQUIRED_FLAGS "-c -std=gnu++14")
set(CMAKE_REQUIRED_QUIET TRUE)

add_custom_target(check-header
    COMMAND ${CMAKE_COMMAND} -E remove CMakeFiles/CMakeError.log
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=ON
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=OFF || true
    COMMAND if [ -f CMakeFiles/CMakeError.log ]\; then cat CMakeFiles/CMakeError.log\; return 1\; else return 0\; fi\;
    WORKING_DIRECTOY ${PROJECT_BINARY_DIR}
    COMMENT "Checking header files"
    USES_TERMINAL
)

# Checks header for standalone compilation
function(_check_header_compilation TARGET)

    get_target_property(SOURCE_FILES ${TARGET} SOURCES)
    get_target_property(SOURCE_DIR ${TARGET} SOURCE_DIR)

    get_directory_property(DEFS DIRECTORY ${SOURCE_DIR} COMPILE_DEFINITIONS)
    foreach(DEF ${DEFS})
        if(${DEF} MATCHES ".*[0-9]\\(.*")
            continue()
        endif()
        list(APPEND DEFS_CLEANED "-D${DEF}")
    endforeach()

    get_target_property(INCLUDE_DIRS ${TARGET} INCLUDE_DIRECTORIES)
    get_target_property(LINK_LIBS ${TARGET} LINK_LIBRARIES)
    foreach(LIB ${LINK_LIBS})
        if(NOT TARGET ${LIB}) # Ignore non-existing targets
            continue()
        endif()
        get_target_property(TARGET_INCLUDE_DIRS ${LIB} INCLUDE_DIRECTORIES)
        if(TARGET_INCLUDE_DIRS)
            list(APPEND INCLUDE_DIRS ${TARGET_INCLUDE_DIRS})
        endif()
    endforeach()
    list(REMOVE_DUPLICATES INCLUDE_DIRS)

    string(REPLACE "${PROJECT_SOURCE_DIR}/" "" DIRECTORY ${SOURCE_DIR})
    message(STATUS "Checking header compilation for ${DIRECTORY} ...")
    include(CheckCXXSourceCompiles)

    set(CMAKE_REQUIRED_INCLUDES ${INCLUDE_DIRS} ${SOURCE_DIR})
    set(CMAKE_REQUIRED_DEFINITIONS ${DEFS_CLEANED})

    foreach(FILE ${SOURCE_FILES})

        if(NOT "${FILE}" MATCHES ".*\\.h") # Check only header files
            continue()
        endif()
        if("${FILE}" MATCHES ".*-impl\\.h") # Ignore *-impl.h files
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
            set(HEADER_COMPILE_ERROR TRUE CACHE INTERNAL "")
            string(REPLACE "${PROJECT_SOURCE_DIR}/" "" FILE_SHORT ${FILE})
            message(STATUS "  Compilation failed for ${FILE_SHORT}")
        endif()
        unset(COMPILES CACHE)

    endforeach()
endfunction()

function(check_header_compilation)
    if(NOT OGS_CHECK_HEADER_COMPILATION)
        return()
    endif()
    set(HEADER_COMPILE_ERROR FALSE CACHE INTERNAL "")

    _check_header_compilation(BaseLib)
    _check_header_compilation(GeoLib)
    _check_header_compilation(MaterialLib)
    _check_header_compilation(MathLib)
    _check_header_compilation(MeshGeoToolsLib)
    _check_header_compilation(MeshLib)
    _check_header_compilation(NumLib)
    _check_header_compilation(ProcessLib)

    if(HEADER_COMPILE_ERROR)
        message(FATAL_ERROR "... header compilation check failed, see CMakeFiles/CMakeError.log for details!")
    endif()
endfunction()
