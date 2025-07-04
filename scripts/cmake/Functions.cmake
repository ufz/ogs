# Returns a list of source files (*.h and *.cpp) in source_files and creates a
# Visual Studio folder. A (relative) subdirectory can be passed as second
# parameter (optional).
macro(GET_SOURCE_FILES source_files)
    if(${ARGC} EQUAL 2)
        set(DIR "${ARGV1}")
    else()
        set(DIR ".")
    endif()

    # Get all files in the directory
    file(GLOB GET_SOURCE_FILES_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
         CONFIGURE_DEPENDS ${DIR}/*.h
    )
    file(GLOB GET_SOURCE_FILES_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
         CONFIGURE_DEPENDS ${DIR}/*.cpp
    )

    set(${source_files}
        ${GET_SOURCE_FILES_HEADERS} ${GET_SOURCE_FILES_TEMPLATES}
        ${GET_SOURCE_FILES_SOURCES}
    )
    list(LENGTH ${source_files} NUM_FILES)
    if(${NUM_FILES} EQUAL 0)
        message(FATAL_ERROR "No source files found in ${DIR}")
    endif()
endmacro()

# Appends a list of source files (*.h and *.cpp) to source_files and creates a
# Visual Studio folder. A (relative) subdirectory can be passed as second
# parameter (optional).
macro(APPEND_SOURCE_FILES source_files)
    if(${ARGC} EQUAL 2)
        set(DIR "${ARGV1}")
    else()
        set(DIR ".")
    endif()

    GET_SOURCE_FILES(TMP_SOURCES "${DIR}")
    set(${source_files} ${${source_files}} ${TMP_SOURCES})
endmacro()

# Adds the include dir containing the autogenerated files to the PUBLIC
# interface of the given target
function(add_autogen_include target)
    get_property(IsMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(IsMultiConfig)
        target_include_directories(
            ${target}
            PUBLIC
                ${CMAKE_CURRENT_BINARY_DIR}/${target}_autogen/include_$<CONFIG>
        )
    else()
        target_include_directories(
            ${target}
            PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${target}_autogen/include
        )
    endif()
endfunction()

# Replacement for add_library() for ogs libraries
function(ogs_add_library targetName)
    set(options STATIC SHARED GENERATE_EXPORT_HEADER)
    cmake_parse_arguments(ogs_add_library "${options}" "" "" ${ARGN})

    foreach(file ${ogs_add_library_UNPARSED_ARGUMENTS})
        get_filename_component(file_path ${file} REALPATH)
        list(APPEND files ${file_path})
    endforeach()

    set(type "") # somehow type may be set to MODULE_LIBRARY from outside, bug?
    if(ogs_add_library_STATIC)
        set(type STATIC)
    elseif(ogs_add_library_SHARED)
        set(type SHARED)
    endif()
    add_library(${targetName} ${type} ${files})
    target_compile_options(
        ${targetName}
        PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Wall -Wextra
                -Wunreachable-code> $<$<CXX_COMPILER_ID:MSVC>:/W3>
    )

    if(BUILD_SHARED_LIBS)
        install(TARGETS ${targetName}
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()

    if(ogs_add_library_GENERATE_EXPORT_HEADER)
        include(GenerateExportHeader)
        generate_export_header(${targetName})
        target_include_directories(
            ${targetName} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
        )
    endif()

    set_target_properties(
        ${targetName} PROPERTIES UNITY_BUILD ${OGS_USE_UNITY_BUILDS}
    )

    if(OGS_INCLUDE_WHAT_YOU_USE)
        if(${CMAKE_CXX_COMPILER_ID} MATCHES ".*Clang")
            set_target_properties(
                ${targetName} PROPERTIES CXX_INCLUDE_WHAT_YOU_USE
                                         include-what-you-use
            )
        else()
            message(
                FATAL_ERROR
                    "OGS_INCLUDE_WHAT_YOU_USE requires the clang compiler!"
            )
        endif()
    endif()

    # Add project root to include directories for cross-library includes
    target_include_directories(${targetName} PRIVATE ${PROJECT_SOURCE_DIR})

    if(MSVC)
        GroupSourcesByFolder(${targetName})
    endif()
endfunction()

set(_OGS_EXECUTABLES "" CACHE INTERNAL "")

# Replacement for ogs_add_executable() for ogs executables
function(ogs_add_executable targetName)
    cmake_parse_arguments(ogs_add_executable "" "" "" ${ARGN})

    foreach(file ${ogs_add_executable_UNPARSED_ARGUMENTS})
        get_filename_component(file_path ${file} REALPATH)
        list(APPEND files ${file_path})
    endforeach()

    if(WIN32)
        # Generate manifest file with long paths enabled. Is added to source
        # files which embeds it into the executable. To work long paths have
        # to be enabled in the registry as well:
        # https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation?tabs=registry#enable-long-paths-in-windows-10-version-1607-and-later
        set(target_mn_filename "${targetName}_longpath.manifest")
        set(mn_file_output "${CMAKE_CURRENT_BINARY_DIR}/${target_mn_filename}")

        set(mn_contents [=[<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<assembly xmlns="urn:schemas-microsoft-com:asm.v1" manifestVersion="1.0">
<application  xmlns="urn:schemas-microsoft-com:asm.v3">
    <windowsSettings xmlns:ws2="http://schemas.microsoft.com/SMI/2016/WindowsSettings">
        <ws2:longPathAware>true</ws2:longPathAware>
    </windowsSettings>
</application>
</assembly>]=])
        file(GENERATE OUTPUT "${mn_file_output}" CONTENT "${mn_contents}")
        list(APPEND files ${mn_file_output})
    endif()

    add_executable(${targetName} ${files})

    target_compile_options(
        ${targetName}
        PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Wall -Wextra
                -Wunreachable-code> $<$<CXX_COMPILER_ID:MSVC>:/W3>
    )

    # Add project root to include directories for cross-library includes
    target_include_directories(${targetName} PRIVATE ${PROJECT_SOURCE_DIR})

    set(ogs_executables "${_OGS_EXECUTABLES}")
    list(APPEND ogs_executables "${targetName}")
    set(_OGS_EXECUTABLES "${ogs_executables}" CACHE INTERNAL "list of exes")
endfunction()

# Parses current directory into a list
function(current_dir_as_list baseDir outList)
    file(RELATIVE_PATH REL_DIR ${PROJECT_SOURCE_DIR}/${baseDir}
         ${CMAKE_CURRENT_LIST_DIR}
    )
    string(REPLACE "/" ";" DIR_LIST ${REL_DIR})
    set(${outList} ${DIR_LIST} PARENT_SCOPE)
endfunction()
