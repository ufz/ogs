# Modified from
# https://github.com/Sbte/BuildExternalProject/commit/ce1a70996aa538aac17a6faf07db487c3a238838
macro(BuildExternalProject_find_package target)

    # Set CMake prefix path so we can look there for the module
    set(_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    mark_as_advanced(_CMAKE_PREFIX_PATH)
    list(APPEND CMAKE_PREFIX_PATH ${build_dir_${target}})

    find_package(${target} MODULE QUIET)
    if(NOT ${target}_FOUND)
        # Look for config version if there was no module
        find_package(
            ${target} CONFIG REQUIRED HINTS ${build_dir_${target}}
            NO_DEFAULT_PATH
        )
    endif()

    # Set CMake prefix path back to what it was
    set(CMAKE_PREFIX_PATH ${_CMAKE_PREFIX_PATH})
    unset(_CMAKE_PREFIX_PATH)
endmacro()

function(BuildExternalProject target)

    message(STATUS "┌─ BuildExternalProject ${target}")
    list(APPEND CMAKE_MESSAGE_INDENT "│    ")

    set(build_dir ${PROJECT_BINARY_DIR}/_ext/${target})
    string(REPLACE ";" " " ARGN_STRING "${ARGN}")

    if(CPM_SOURCE_CACHE)
        cmake_path(
            IS_PREFIX PROJECT_BINARY_DIR "${CPM_SOURCE_CACHE}" _is_inside_build
        )
        if(NOT _is_inside_build)
            string(SHA256 _hash "${CMAKE_GENERATOR};${ARGN_STRING}")
            set(build_dir "${CPM_SOURCE_CACHE}/_ext/${target}/${_hash}")
            file(LOCK ${build_dir}/cmake.lock)
        endif()
    endif()

    message(STATUS "Building ${target} in ${build_dir}")
    set(build_dir_${target} "${build_dir}" CACHE INTERNAL "")

    file(MAKE_DIRECTORY ${build_dir})

    set(CMAKE_LIST_CONTENT
        "
    cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})
    project(externalproject_${target})

    include(ExternalProject)
    ExternalProject_add(${target}
      PREFIX ${build_dir}
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      LOG_DOWNLOAD YES
      LOG_UPDATE YES
      LOG_PATCH YES
      LOG_CONFIGURE YES
      LOG_BUILD YES
      LOG_TEST YES
      LOG_INSTALL YES
      LOG_MERGED_STDOUTERR YES
      LOG_DIR logs
      ${ARGN_STRING}
      )

    add_custom_target(build_${target})
    add_dependencies(build_${target} ${target})
    "
    )

    if(EXISTS ${build_dir}/CMakeLists.txt)
        file(SHA256 ${build_dir}/CMakeLists.txt file_sha)
        string(SHA256 new_sha "${CMAKE_LIST_CONTENT}")

        if(NOT file_sha STREQUAL new_sha)
            file(WRITE ${build_dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
            BuildExternalProject_configure(${build_dir})
        endif()
    else()
        file(WRITE ${build_dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
        BuildExternalProject_configure(${build_dir})
    endif()

    message(
        STATUS
            "Finished building ${target}. Logs in ${build_dir}/src/${target}-stamp"
    )
    if(EXISTS ${build_dir}/cmake.lock)
        file(LOCK ${build_dir}/cmake.lock RELEASE)
    endif()
    list(POP_BACK CMAKE_MESSAGE_INDENT)
    message(STATUS "└─ End BuildExternalProject ${_target}")
endfunction()

function(BuildExternalProject_configure build_dir)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result WORKING_DIRECTORY ${build_dir}
    )

    if(result)
        message(FATAL_ERROR "CMake step for external project failed: ${result}")
    else()
        BuildExternalProject_build(${build_dir})
    endif()
endfunction()

function(BuildExternalProject_build build_dir)
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
        RESULT_VARIABLE result WORKING_DIRECTORY ${build_dir}
    )

    if(result)
        message(FATAL_ERROR "Build step for external project failed: ${result}")
    endif()
endfunction()
