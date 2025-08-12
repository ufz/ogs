if(LINUX)
    add_executable(ogs_test_dlopen dlopen.cpp)

    # Adds ogs_test_dlopen tests for each shared library in the passed dir.
    # Any surplus arguments (${ARGN}) are added as commandline arguments to the test command.
    function(add_dlopen_tests dir)
        get_directory_property(_targets DIRECTORY "${dir}" BUILDSYSTEM_TARGETS)

        foreach(target IN LISTS _targets)
            get_target_property(_type ${target} TYPE)
            if(_type STREQUAL "SHARED_LIBRARY")
                add_test(
                    NAME "dlopen_${target}"
                    COMMAND $<TARGET_FILE:ogs_test_dlopen> ${ARGN} $<TARGET_FILE:${target}>)
                set_tests_properties("dlopen_${target}" PROPERTIES DEPENDS ogs_test_dlopen)
            endif()
        endforeach()
    endfunction()

    # Invokes the passed callback for each subdirectory of the passed directory (and for dir itself).
    # Only subdirectories of the ${PROJECT_SOURCE_DIR} are considered.
    # Any surplus arguments (${ARGN}) are passed on to the callback.
    function(for_each_subdir dir callback)
        if(COMMAND ${callback})
            cmake_language(CALL "${callback}" "${dir}" ${ARGN})
        else()
            message(FATAL_ERROR "Callback ${callback} is not a valid function")
        endif()

        get_directory_property(_subdirs DIRECTORY "${dir}" SUBDIRECTORIES)

        foreach(subdir IN LISTS _subdirs)
            string(FIND "${subdir}" "${PROJECT_SOURCE_DIR}" _index)
            if(NOT _index EQUAL 0)
                continue()
            endif()

            for_each_subdir("${subdir}" "${callback}" ${ARGN})
        endforeach()
    endfunction()
endif()
