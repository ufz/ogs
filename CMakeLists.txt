# ---- OGS-6 Project ----
cmake_minimum_required(VERSION 3.31)
cmake_policy(SET CMP0155 OLD) # Disable module import scanning.

if(APPLE AND NOT (DEFINED CMAKE_OSX_DEPLOYMENT_TARGET
                  OR DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
)
    # Setting the deployment target to current macOS major version
    cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
    string(REGEX MATCH "^[0-9]+" _os_release_major "${_os_release}")
    set(CMAKE_OSX_DEPLOYMENT_TARGET ${_os_release_major}.0
        CACHE STRING "Minimum macOS operating system version for deployment"
    )
endif()

project(OGS-6)

include(scripts/cmake/CMakeSetup.cmake)
include(ParseCMakeArgs)
include(Options)

include(scripts/cmake/UseCompilerCache.cmake)
if(NOT OGS_DISABLE_COMPILER_CACHE)
    useCompilerCache()
endif()
# Suppress dev-warning in CPMFindPackage(NAME Boost)
set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
include(scripts/cmake/CPM.cmake)
if(DEFINED ENV{CI})
    message(STATUS "CMake arguments: ${CMAKE_ARGS}")
endif()
include(CTest)

# ---- CMake includes ----
include(Versions)
include(GitSetup)
include(GuixSetup)
include(PythonSetup)
include(ProcessesSetup)
include(ProjectSetup)
include(Functions)
include(CompilerSetup)
include(JobPools)
include(Find)
include(DependenciesExternalProject)
include(Dependencies)
include(DocumentationSetup)
if(OGS_BUILD_TESTING)
    include(test/Test)
endif()
include(Coverage)
include(CppCheck)
include(WebSetup)

set(CMAKE_LINK_LIBRARIES_ONLY_TARGETS ON)

include(packaging/Pack)

# ---- Subdirectories ----
add_subdirectory(Applications)
add_subdirectory(BaseLib)
add_subdirectory(GeoLib)
add_subdirectory(InfoLib)
add_subdirectory(MathLib)
add_subdirectory(MeshLib)
add_subdirectory(MeshToolsLib)
add_subdirectory(MeshGeoToolsLib)
add_subdirectory(NumLib)

if(_build_chemistry_lib)
    add_subdirectory(ChemistryLib)
endif()

if(OGS_BUILD_CLI OR OGS_BUILD_UTILS OR OGS_BUILD_TESTING)
    if(TARGET OgsXdmf OR Xdmf_FOUND)
        add_subdirectory(Tests/xdmfdiff)
    endif()
    add_subdirectory(ParameterLib)
    add_subdirectory(MaterialLib)
    add_subdirectory(ProcessLib)
endif()
if(OGS_BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
    add_subdirectory(Tests)
endif()

include(UnityBuildSettings)
include(GenerateEnvrcFile)
# Add an alias target to all_verify_interface_header_sets. Individual libs can
# be checked, e.g.: MathLib_verify_interface_header_sets
add_custom_target(check-headers)
add_dependencies(check-headers all_verify_interface_header_sets)
include(Features)
printEnabledProcesses()
printMKLUsage()

# Check wheel binaries
list(SORT _OGS_EXECUTABLES)
# SWMMConverter is win-only
list(REMOVE_ITEM _OGS_EXECUTABLES testrunner SWMMConverter)

if(OGS_BUILD_WHEEL)
    execute_process(
        COMMAND ${Python_EXECUTABLE} scripts/python/compare_binaries.py
                "${_OGS_EXECUTABLES}" ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
        RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
        message(
            FATAL_ERROR
                "Binary list comparison failed. Please update binaries_list in provide_ogs_cli_tools_via_wheel.py!"
        )
    endif()
endif()

if(TARGET ogs_test_dlopen)
    # Has to come last in this top-level CMakeLists.txt to add tests for all
    # shared libraries. We pass ${Python_LIBRARIES} as additional lib(s) to
    # load, because many OGS libs (e.g. each process, because of a cyclic
    # dependency between the generic process lib and (Python) BCs) depend on
    # them and they are not explicitly linked against.
    for_each_subdir(
        "${PROJECT_SOURCE_DIR}" add_dlopen_tests ${Python_LIBRARIES}
    )
endif()
