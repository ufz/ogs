if(USE_CONAN)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    include_directories(SYSTEM ${CONAN_INCLUDE_DIRS_VTK})
    return()
endif()
include(ThirdPartyLibVersions)
include(ExternalProject)

if(NOT DEFINED VTK_DIR AND DEFINED ENV{VTK_DIR})
    set(VTK_DIR $ENV{VTK_DIR})
endif()

set(VTK_LIBRARIES ${VTK_MODULES} CACHE STRING "" FORCE)
if(OGS_BUILD_GUI)
    # Replace vtknetcdf with vtkNetCDF vtkNetCDF_cxx
    list(REMOVE_ITEM VTK_LIBRARIES vtknetcdf)
    list(APPEND VTK_LIBRARIES vtkNetCDF vtkNetCDF_cxx)
endif()

if(OGS_LIB_VTK STREQUAL "System")
    find_package(VTK COMPONENTS ${VTK_MODULES} NO_MODULE REQUIRED)
    if(NOT VTK_FOUND)
        message(FATAL_ERROR "Aborting CMake because system VTK was not found!")
    endif()
elseif(OGS_LIB_VTK STREQUAL "Default" OR DEFINED VTK_DIR)
    find_package(VTK COMPONENTS ${VTK_MODULES} NO_MODULE QUIET)
endif()

if(VTK_FOUND)
    if(VTK_VERSION VERSION_LESS ${OGS_VTK_VERSION})
        message(FATAL_ERROR "Aborting CMake because VTK ${VTK_VERSION} is too old! (required: ${OGS_VTK_VERSION})")
    endif()
    message(STATUS "Using VTK ${VTK_VERSION} in ${VTK_DIR}")
    foreach(DIR ${VTK_INCLUDE_DIRS})
            if("${DIR}" MATCHES ".*vtknetcdf.*")
                include_directories(SYSTEM ${DIR}/../cxx ${DIR}/include)
            elseif("${DIR}" MATCHES ".*vtk.*")
                include_directories(SYSTEM ${DIR}/vtknetcdf/include)
            endif()
        endforeach()
        include_directories(SYSTEM ${VTK_DIR}/../ThirdParty/netcdf/vtknetcdf/cxx)
    return()
endif()
set(VTK_DIR ${CMAKE_BINARY_DIR}/External/vtk/src/vtk-build CACHE PATH "" FORCE)

message(STATUS "Building VTK as an external project in the build directory")

if(WIN32 AND NOT CYGWIN)
    set(VTK_MAKE_COMMAND
        msbuild /p:Configuration=Release /m:${NUM_PROCESSORS} VTK.sln &&
        msbuild /p:Configuration=Debug /m:${NUM_PROCESSORS} /m VTK.sln)
else()
    if($ENV{CI})
        set(VTK_MAKE_COMMAND make)
    else()
        set(VTK_MAKE_COMMAND make -j ${NUM_PROCESSORS})
    endif()
endif()

# Enable just the modules we selected
set(VTK_CMAKE_ARGS -DVTK_Group_StandAlone:bool=off -DVTK_Group_Rendering:bool=off)
foreach(arg ${VTK_MODULES})
    list(APPEND VTK_CMAKE_ARGS -DModule_${arg}:bool=on)
endforeach()

ExternalProject_Add(vtk
    PREFIX ${CMAKE_BINARY_DIR}/External/vtk
    URL ${OGS_VTK_URL}
    CMAKE_ARGS -Wno-dev
    CMAKE_CACHE_ARGS ${VTK_CMAKE_ARGS}
    BUILD_COMMAND ${VTK_MAKE_COMMAND}
    INSTALL_COMMAND ""
)

if(NOT VTK_FOUND)
    # Rerun cmake in initial build
    add_custom_target(VtkRescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} DEPENDS vtk)
else()
    add_custom_target(VtkRescan) # dummy target for caching
endif()
