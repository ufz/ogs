# This file initializes the required submodules
set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "" CACHE STRING
    "User given submodules which should be checked out by CMake.")
if(NOT OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT)
    set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "")
endif()
set(REQUIRED_SUBMODULES
    ThirdParty/autocheck
    ThirdParty/cmake-modules
    ThirdParty/vtkdiff
    ThirdParty/tclap
    ThirdParty/tetgen
    ${OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT}
)
if(OGS_BUILD_TESTS)
    list(APPEND REQUIRED_SUBMODULES Tests/Data)
endif()
if(OGS_BUILD_GUI)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/vtkGUISupportQt)
endif()
if(OGS_BUILD_METIS)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/metis)
endif()
if(OGS_BUILD_SWMM)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/SwmmInterface)
endif()

# Sync submodules, which is required when a submodule changed its URL
if(OGS_SYNC_SUBMODULES)
    execute_process(
        COMMAND ${GIT_TOOL_PATH} submodule sync
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_QUIET
    )
endif()
foreach(SUBMODULE ${REQUIRED_SUBMODULES})
    execute_process(
        COMMAND ${GIT_TOOL_PATH} submodule status ${SUBMODULE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE SUBMODULE_STATE
    )
    string(REGEX MATCH "^\\-" UNINITIALIZED ${SUBMODULE_STATE})
    string(REGEX MATCH "^\\+" MISMATCH ${SUBMODULE_STATE})

    if(UNINITIALIZED)
        message(STATUS "Initializing submodule ${SUBMODULE}")
        execute_process(
            COMMAND ${GIT_TOOL_PATH}
                submodule update --init --recursive ${SUBMODULE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    elseif(MISMATCH)
        message(STATUS "Updating submodule ${SUBMODULE}")
        execute_process(
            COMMAND ${GIT_TOOL_PATH}
                submodule update --recursive ${SUBMODULE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()
endforeach()
