# Supply include directories and compiler flags
get_directory_property(INCLUDE_DIRS INCLUDE_DIRECTORIES)
set(CMAKE_REQUIRED_FLAGS "-c")

add_custom_target(check-header
    COMMAND ${CMAKE_COMMAND} -E remove -f CMakeFiles/CMakeError.log
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=ON
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=OFF || true
    COMMAND if [ -f CMakeFiles/CMakeError.log ]\; then cat CMakeFiles/CMakeError.log\; return 1\; else return 0\; fi\;
    WORKING_DIRECTOY ${PROJECT_BINARY_DIR}
    COMMENT "Checking header files"
    USES_TERMINAL
)

# Checks header for standalone compilation
function(_check_header_compilation TARGET)

    if(NOT TARGET ${TARGET})
        return()
    endif()

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
        # Ignore non-existing targets or interface libs
        if(NOT TARGET ${LIB})
            continue()
        endif()
        get_target_property(LIB_TYPE ${LIB} TYPE)
        if(LIB_TYPE STREQUAL "INTERFACE_LIBRARY")
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
    # HACK, maybe add Gui Widgets Xml XmlPatterns as well
    if(OGS_BUILD_GUI)
        set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES}
            ${Qt5Core_INCLUDE_DIRS}
            ${Qt5Gui_INCLUDE_DIRS}
            ${Qt5Widgets_INCLUDE_DIRS}
        )
    endif()
    set(CMAKE_REQUIRED_DEFINITIONS ${DEFS_CLEANED})

    foreach(FILE ${SOURCE_FILES})

        if(NOT "${FILE}" MATCHES ".*\\.h") # Check only header files
            continue()
        endif()
        if("${FILE}" MATCHES ".*-impl\\.h") # Ignore *-impl.h files
            continue()
        endif()
        if("${FILE}" MATCHES ".*Dialog\\.h") # Ignore Qt Dialog files
            continue()
        endif()
        if("${FILE}" MATCHES ".*Widget\\.h") # Ignore Qt Widget files
            continue()
        endif()
        if("${FILE}" MATCHES ".*Window\\.h") # Ignore Qt Window files
            continue()
        endif()

        string(REPLACE "${PROJECT_SOURCE_DIR}/" "" TEST_NAME ${FILE})
        string(REPLACE "." "_" TEST_NAME ${TEST_NAME})
        string(REPLACE "/" "_" TEST_NAME ${TEST_NAME})
        check_cxx_source_compiles(
            "
            #include \"${FILE}\"
            int main() { return 0; }
            "
            ${TEST_NAME}_COMPILES
        )

        if(NOT ${TEST_NAME}_COMPILES)
            set(HEADER_COMPILE_ERROR TRUE CACHE INTERNAL "")
            message(STATUS "  Compilation failed for ${FILE}")
        endif()
        unset(${TEST_NAME}_COMPILES CACHE)

        unset(TEST_NAME)
    endforeach()
endfunction()

function(check_header_compilation)
    if(NOT OGS_CHECK_HEADER_COMPILATION)
        return()
    endif()
    set(HEADER_COMPILE_ERROR FALSE CACHE INTERNAL "")

    _check_header_compilation(BaseLib)
    _check_header_compilation(ChemistryLib)
    _check_header_compilation(GeoLib)
    foreach(lib Git CMake Test)
        _check_header_compilation(${lib}InfoLib)
    endforeach(lib)
    _check_header_compilation(MaterialLib)
    _check_header_compilation(MathLib)
    _check_header_compilation(MeshGeoToolsLib)
    _check_header_compilation(MeshLib)
    _check_header_compilation(NumLib)
    _check_header_compilation(ParameterLib)
    _check_header_compilation(ProcessLib)
    _check_header_compilation(ApplicationsLib)
    _check_header_compilation(ApplicationsFileIO)
    _check_header_compilation(DataHolderLib)
    if(OGS_BUILD_GUI)
        _check_header_compilation(QtBase)
        _check_header_compilation(QtDataView)
        _check_header_compilation(QtDiagramView)
        # _check_header_compilation(QtStratView) # all fail
        _check_header_compilation(VtkVis)
    endif()

    if(HEADER_COMPILE_ERROR)
        message(FATAL_ERROR "... header compilation check failed, see CMakeFiles/CMakeError.log for details!")
    endif()
endfunction()
