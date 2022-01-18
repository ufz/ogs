# cmake-lint: disable=R0912,R0915,C0103

# Supply include directories and compiler flags
get_directory_property(INCLUDE_DIRS INCLUDE_DIRECTORIES)
set(CMAKE_REQUIRED_FLAGS "-c -fPIC")

add_custom_target(
    check-header
    COMMAND ${CMAKE_COMMAND} -E remove -f CMakeFiles/CMakeError.log
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=ON
    COMMAND ${CMAKE_COMMAND} . -DOGS_CHECK_HEADER_COMPILATION=OFF || true
    COMMAND
        if [ -f CMakeFiles/CMakeError.log ]\; then cat
        CMakeFiles/CMakeError.log\; return 1\; else return 0\; fi\;
        WORKING_DIRECTOY ${PROJECT_BINARY_DIR}
    COMMENT "Checking header files"
    USES_TERMINAL
)

# Checks header for standalone compilation
function(_check_header_compilation target)
    # cmake-lint: disable=R0915
    if(NOT TARGET ${target})
        return()
    endif()

    get_target_property(SOURCE_FILES ${target} SOURCES)
    get_target_property(SOURCE_DIR ${target} SOURCE_DIR)

    get_directory_property(DEFS DIRECTORY ${SOURCE_DIR} COMPILE_DEFINITIONS)
    foreach(def ${DEFS})
        if(${def} MATCHES ".*[0-9]\\(.*")
            continue()
        endif()
        list(APPEND DEFS_CLEANED "-D${def}")
    endforeach()

    get_target_property(INCLUDE_DIRS ${target} INCLUDE_DIRECTORIES)
    get_target_property(
        INTERFACE_INCLUDE_DIRS ${target} INTERFACE_INCLUDE_DIRECTORIES
    )
    if(INTERFACE_INCLUDE_DIRS)
        list(APPEND INCLUDE_DIRS ${INTERFACE_INCLUDE_DIRS})
    endif()
    get_target_property(LINK_LIBS ${target} LINK_LIBRARIES)
    # Transitive dependencies are not resolved
    foreach(lib ${LINK_LIBS} spdlog::spdlog Boost::boost Eigen3::Eigen
                nlohmann_json::nlohmann_json
    )
        # Ignore non-existing targets or interface libs
        if(NOT TARGET ${lib})
            continue()
        endif()
        get_target_property(LIB_TYPE ${lib} TYPE)
        if(LIB_TYPE STREQUAL "INTERFACE_LIBRARY")
            get_target_property(
                TARGET_INCLUDE_DIRS ${lib} INTERFACE_INCLUDE_DIRECTORIES
            )
        else()
            get_target_property(TARGET_INCLUDE_DIRS ${lib} INCLUDE_DIRECTORIES)
        endif()
        if(TARGET_INCLUDE_DIRS)
            if("${TARGET_INCLUDE_DIRS}" MATCHES ".*<BUILD_INTERFACE:([^>]*)")
                list(APPEND INCLUDE_DIRS ${CMAKE_MATCH_1})
            else()
                list(APPEND INCLUDE_DIRS ${TARGET_INCLUDE_DIRS})
            endif()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES INCLUDE_DIRS)

    string(REPLACE "${PROJECT_SOURCE_DIR}/" "" DIRECTORY ${SOURCE_DIR})
    message(STATUS "Checking header compilation for ${DIRECTORY} ...")
    include(CheckCXXSourceCompiles)

    # cmake-lint: disable=C0103
    set(CMAKE_REQUIRED_INCLUDES ${INCLUDE_DIRS} ${SOURCE_DIR})
    # HACK, maybe add Gui Widgets Xml XmlPatterns as well
    if(OGS_BUILD_GUI)
        set(CMAKE_REQUIRED_INCLUDES
            ${CMAKE_REQUIRED_INCLUDES} ${Qt5Core_INCLUDE_DIRS}
            ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS}
        )
    endif()
    get_target_property(COMPILE_DEFINITIONS ${target} _target_defs)
    set(CMAKE_REQUIRED_DEFINITIONS ${_target_defs} ${DEFS_CLEANED})

    foreach(file ${SOURCE_FILES})

        if(NOT "${file}" MATCHES ".*\\.h") # Check only header files
            continue()
        endif()
        if("${file}" MATCHES ".*-impl\\.h") # Ignore *-impl.h files
            continue()
        endif()
        if("${file}" MATCHES ".*Dialog\\.h") # Ignore Qt Dialog files
            continue()
        endif()
        if("${file}" MATCHES ".*Widget\\.h") # Ignore Qt Widget files
            continue()
        endif()
        if("${file}" MATCHES ".*Window\\.h") # Ignore Qt Window files
            continue()
        endif()
        if("${file}" MATCHES "ui_.*\\.h") # Ignore Qt-generated ui files
            continue()
        endif()

        string(REPLACE "${PROJECT_SOURCE_DIR}/" "" TEST_NAME ${file})
        string(REPLACE "." "_" TEST_NAME ${TEST_NAME})
        string(REPLACE "/" "_" TEST_NAME ${TEST_NAME})
        check_cxx_source_compiles(
            "
            #include \"${file}\"
            int main() { return 0; }
            "
            ${TEST_NAME}_COMPILES
        )

        if(NOT ${TEST_NAME}_COMPILES)
            set(_HEADER_COMPILE_ERROR TRUE CACHE INTERNAL "")
            message(STATUS "  Compilation failed for ${file}")
        endif()
        unset(${TEST_NAME}_COMPILES CACHE)

        unset(TEST_NAME)
    endforeach()
endfunction()

# Check header compilation in
function(check_header_compilation)
    if(NOT OGS_CHECK_HEADER_COMPILATION)
        return()
    endif()
    set(_HEADER_COMPILE_ERROR FALSE CACHE INTERNAL "")

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

    if(_HEADER_COMPILE_ERROR)
        message(
            FATAL_ERROR
                "... header compilation check failed, see CMakeFiles/CMakeError.log for details!"
        )
    endif()
endfunction()
