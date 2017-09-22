macro(InstallDependencies TARGET)
    set(EXCLUDE_SYSTEM 1)
    get_target_property(EXE_DIR ${TARGET} RUNTIME_OUTPUT_DIRECTORY)
    set(TARGET_EXE ${EXE_DIR}/${TARGET}${CMAKE_EXECUTABLE_SUFFIX})

    if(EXISTS ${TARGET_EXE})
        # Run CMake again after target was built to collect dependencies
        add_custom_command(TARGET ${TARGET} POST_BUILD
            COMMAND ${CMAKE_COMMAND} . -DPRE_INSTALL_RUN=ON
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )

        include(GetPrerequisites)
        # arg3: exclude system, arg4: recursive
        if (VTK_BUILD_SHARED_LIBS)
            list(APPEND dirs ${vtkIOXML_RUNTIME_LIBRARY_DIRS})
        endif()
        list(APPEND dirs "/usr/local/lib")
        get_prerequisites(${TARGET_EXE} TARGET_DEPENDENCIES ${EXCLUDE_SYSTEM} 1 "" ${dirs})
        if(PRE_INSTALL_RUN)
            message("-- Dependencies of target ${TARGET}:")
        endif()
        foreach(DEPENDENCY ${TARGET_DEPENDENCIES})
            if(NOT ${DEPENDENCY} MATCHES "@loader_path")
                gp_resolve_item("${TARGET_EXE}" "${DEPENDENCY}" "" "" DEPENDENCY_PATH)
                get_filename_component(RESOLVED_DEPENDENCY_PATH "${DEPENDENCY_PATH}" REALPATH)
                string(TOLOWER ${DEPENDENCY} DEPENDENCY_LOWER)
                set(DEPENDENCY_PATHS ${DEPENDENCY_PATHS} ${RESOLVED_DEPENDENCY_PATH})
                if(PRE_INSTALL_RUN)
                    message("     ${RESOLVED_DEPENDENCY_PATH}")
                endif()
            endif()
        endforeach()
        if(PRE_INSTALL_RUN)
            message("")
        endif()
        install(FILES ${DEPENDENCY_PATHS} DESTINATION bin)
    endif()

endmacro()

