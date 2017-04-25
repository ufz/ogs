if(USE_CONAN)
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${CONAN_INCLUDE_DIRS_EIGEN3})
    return()
endif()

if(OGS_LIB_EIGEN STREQUAL "System")
    find_package(Eigen3 3.2.9 REQUIRED)
    if(NOT EIGEN3_FOUND)
        message(FATAL_ERROR "Aborting CMake because system Eigen was not found!")
    endif()
elseif(OGS_LIB_EIGEN STREQUAL "Default")
    find_package(Eigen3 3.2.9)
endif()

# First check for system Eigen
if(NOT EIGEN3_INCLUDE_DIR)
    if(EIGEN3_FOUND)
        set(EIGEN3_FOUND TRUE CACHE BOOL "Was Eigen found?" FORCE)
        set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR}" CACHE STRING "Eigen include dir" FORCE)
        return()
    else()
        set(EIGEN3_INCLUDE_DIR "")
    endif()
endif()

if(EIGEN3_FOUND)
    return()
endif()

include(ThirdPartyLibVersions)
include(ExternalProject)
ExternalProject_Add(Eigen
    PREFIX ${PROJECT_BINARY_DIR}/External/eigen
    URL ${OGS_EIGEN_URL}
    URL_MD5 ${OGS_EIGEN_MD5}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property( Eigen source_dir )

if(NOT EIGEN3_INCLUDE_DIR)
    set( EIGEN3_INCLUDE_DIR ${source_dir} CACHE INTERNAL "Eigen include dir" FORCE)
    message(STATUS "Downloading Eigen automatically.")
endif()
