message(STATUS "********* Conan FindBoost wrapper **********")
message("COMPONENTS TO SEARCH: ${Boost_FIND_COMPONENTS}")

set(BOOST_ROOT ${CONAN_BOOST_ROOT})
set(BOOST_INCLUDEDIR ${CONAN_INCLUDE_DIRS_BOOST})
set(Boost_LIBRARY_DIR ${CONAN_LIB_DIRS_BOOST})
set(BOOST_LIBRARYDIR ${CONAN_LIB_DIRS_BOOST})
set(Boost_NO_SYSTEM_PATHS ON)
set(Boost_NO_BOOST_CMAKE ON)

# READ conaninfo and detect HEADER ONLY
file(READ ${CONAN_BOOST_ROOT}/conaninfo.txt CONANINFO_FILE)
if(WIN32)
    # Appends "g"
    if(CONANINFO_FILE MATCHES "build_type=Debug")
        set(Boost_USE_DEBUG_RUNTIME ON)
    else()
        set(Boost_USE_DEBUG_RUNTIME OFF)
    endif()

    # Appends "s"
    if(CONANINFO_FILE MATCHES "compiler.runtime=MT" OR CONANINFO_FILE MATCHES "compiler.runtime=MTd")
        set(Boost_USE_STATIC_RUNTIME ON)
    else()
        set(Boost_USE_STATIC_RUNTIME OFF)
    endif()

    message("DEBUG RUNTIME: ${Boost_USE_DEBUG_RUNTIME}")
    message("STATIC RUNTIME: ${Boost_USE_STATIC_RUNTIME}")

    # The space is important, so it doesn't match the flag for zlib:shared=False
    if(CONANINFO_FILE MATCHES " shared=False")
        set(Boost_USE_STATIC_LIBS ON) # Removed in the original file
    else()
        set(Boost_USE_STATIC_LIBS OFF)
    endif()
endif()

include(${CMAKE_SOURCE_DIR}/scripts/cmake/conan/OriginalFindBoost_3_4_3.cmake)
