if(NOT IS_GIT_REPO)
    return()
endif()

# This file initializes the required submodules
set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "" CACHE STRING
    "User given submodules which should be checked out by CMake.")
if(NOT OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT)
    set(OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT "")
endif()
set(REQUIRED_SUBMODULES
    ThirdParty/autocheck
    ThirdParty/cmake-modules
    ThirdParty/exprtk
    ThirdParty/googletest
    ThirdParty/iphreeqc/src
    ThirdParty/json-cmake
    ThirdParty/spdlog
    ThirdParty/tclap
    ThirdParty/tetgen
    ${OGS_ADDITIONAL_SUBMODULES_TO_CHECKOUT}
)
if(BUILD_TESTING)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/vtkdiff)
endif()
if(OGS_BUILD_UTILS)
    # Required by the partmesh tool, which is build with utils only.
    list(APPEND REQUIRED_SUBMODULES ThirdParty/metis)
endif()
if(OGS_BUILD_SWMM)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/SwmmInterface)
endif()
if(OGS_USE_PYTHON)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/pybind11)
endif()
if (OGS_USE_MFRONT)
    list(APPEND REQUIRED_SUBMODULES ThirdParty/MGIS)
endif()

# Sync submodules, which is required when a submodule changed its URL
if(OGS_SYNC_SUBMODULES)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule sync
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_QUIET
    )
endif()
foreach(SUBMODULE ${REQUIRED_SUBMODULES})
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule status ${SUBMODULE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE SUBMODULE_STATE
    )
    string(REGEX MATCH "^\\-" UNINITIALIZED ${SUBMODULE_STATE})
    string(REGEX MATCH "^\\+" MISMATCH ${SUBMODULE_STATE})

    if(IS_CI)
        # Always set submodule to the given state
        execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init --force
                --recursive ${DEPTH} ${SUBMODULE}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            RESULT_VARIABLE RESULT
        )
    else()
        set(RESULT "")
        if(UNINITIALIZED)
            message(STATUS "Initializing submodule ${SUBMODULE}")
            execute_process(
                COMMAND ${GIT_EXECUTABLE} submodule update --init
                    --recursive ${DEPTH} ${SUBMODULE}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE RESULT
            )

        elseif(MISMATCH)
            message(STATUS "Updating submodule ${SUBMODULE}")
            execute_process(
                COMMAND ${GIT_EXECUTABLE} submodule update
                    --recursive ${SUBMODULE}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE RESULT
            )
        endif()
    endif()

    if((NOT ${RESULT} STREQUAL "") AND (NOT ${RESULT} STREQUAL "0"))
        message(FATAL_ERROR "Error in submodule setup; return value: ${RESULT}")
    endif()
endforeach()
