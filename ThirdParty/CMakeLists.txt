# Disable all warnings
if(COMPILER_IS_CLANG OR COMPILER_IS_GCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
endif()

# Add logog subdirectory and group its targets in a Visual Studio folder
add_subdirectory(logog)
if(BUILD_TESTING)
    set_property(TARGET Continuous PROPERTY FOLDER "logog")
endif()
set_property(TARGET logog PROPERTY FOLDER "logog")
set_property(TARGET test-logog PROPERTY FOLDER "logog")
set_target_properties(test-logog PROPERTIES EXCLUDE_FROM_ALL TRUE)
if(MSVC)
    set_target_properties(logog PROPERTIES COMPILE_FLAGS /W0)
endif()

add_subdirectory(gtest)

# This is a workaround for Travis builds.
option(OGS_USE_OPTIONAL_SUBMODULES "Option for enabling optional submodules" OFF)

# VtkFbxConverter
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VtkFbxConverter/CMakeLists.txt
    AND OGS_USE_OPTIONAL_SUBMODULES)
    add_subdirectory(VtkFbxConverter)
endif()

# tetgen
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tetgen/CMakeLists.txt)
    add_subdirectory(tetgen)
    install (PROGRAMS $<TARGET_FILE:tetgen> DESTINATION bin COMPONENT ogs_extras)
endif()

# vtkdiff
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/vtkdiff/CMakeLists.txt AND OGS_BUILD_TESTS)
    add_subdirectory(vtkdiff)
    install (PROGRAMS $<TARGET_FILE:vtkdiff> DESTINATION bin COMPONENT ogs_extras)
endif()

if(OGS_BUILD_SWMM)
    add_subdirectory(SWMMInterface)
    set_target_properties(SWMM SwmmInterface PROPERTIES COMPILE_FLAGS /W0)
endif()

if(OGS_BUILD_METIS)
    include(${PROJECT_SOURCE_DIR}/scripts/cmake/MetisSetup.cmake)
endif()
