if(USE_CONAN)
    SET(BOOST_HEADER_ONLY TRUE)
    find_package(Boost REQUIRED)
    include_directories(SYSTEM ${CONAN_INCLUDE_DIRS_BOOST})
    link_directories(${Boost_LIBRARY_DIR})
    return()
endif()
include(ThirdPartyLibVersions)
include(ExternalProject)

if(Boost_FOUND)
    return()
endif()

if(NOT DEFINED BOOST_ROOT AND DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "")
endif()

# First check for system boost
if(NOT Boost_INCLUDE_DIRS)
    if(APPLE)
        set(BOOST_ROOT $ENV{HOMEBREW_ROOT})
    endif()
    if(WIN32 AND COMPILER_IS_GCC)
        set(BOOST_INCLUDEDIR "$ENV{CMAKE_LIBRARY_SEARCH_PATH}/include/boost*")
    endif()
    if(OGS_LIB_BOOST STREQUAL "System")
        find_package(Boost ${OGS_BOOST_VERSION} REQUIRED)
        if(NOT Boost_FOUND)
            message(FATAL_ERROR "Aborting CMake because system Boost was not found!")
        endif()
    elseif(OGS_LIB_BOOST STREQUAL "Default")
        find_package(Boost ${OGS_BOOST_VERSION})
    endif()
    if(Boost_FOUND)
        set(Boost_FOUND TRUE CACHE BOOL "Was Boost found?" FORCE)
        set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}" CACHE STRING "Boost include dirs" FORCE)
        return()
    else()
        set(Boost_INCLUDE_DIRS "")
    endif()
endif()

ExternalProject_Add(Boost
    PREFIX ${PROJECT_BINARY_DIR}/External/boost
    URL ${OGS_BOOST_URL}
    URL_MD5 ${OGS_BOOST_MD5}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property( Boost source_dir )

if(NOT Boost_INCLUDE_DIRS)
    set( Boost_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Boost include directories")
    message(STATUS "Downloading Boost ${OGS_BOOST_VERSION} automatically.")
endif()
